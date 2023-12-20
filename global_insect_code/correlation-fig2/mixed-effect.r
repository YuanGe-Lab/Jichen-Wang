#code for mixed effects model
library(data.table)
li <- c("weighted_unifrac_distance_matrix")
i <- 1
bray <- fread(paste0("../../core-metrics-results/",li[i],"/distance-matrix.tsv"),
              header = T)
bray <- as.data.frame(bray)
rownames(bray) <- bray$V1
bray <- bray[,-1]
sum(rownames(bray) == colnames(bray))

meta <- read.delim("../../meta.txt")
meta <- meta[,c(1,16,30,31)]
meta <- meta[meta$order != "Primates",]
rownames(meta) <- meta$X
sum(rownames(meta) == rownames(bray))
order <- unique(meta$order)

di <- c("ML_wd.txt","NJ_wd.txt")
j <- 1
phydis <- fread(paste0("../../new_15928samples/tree/",di[j]),
                header = T)
phydis <- as.data.frame(phydis)
rownames(phydis) <- colnames(phydis)
sum(colnames(phydis) == colnames(bray))

k <- 3
gl <- c("site_dis.txt","dis_2.txt","dis_whole.txt")
geocl <- fread(paste0("../../map_geodistance/distance2_whole/",gl[k]),
               header = T)
geocl <- as.data.frame(geocl)
rownames(geocl) <- geocl$V1
geocl <- geocl[,-1]
sum(rownames(geocl) == rownames(bray))

k <- 1
geo <- fread(paste0("../../map_geodistance/distance2_whole/",gl[k]),
             header = T)
geo <- as.data.frame(geo)
rownames(geo) <- geo$V1
geo <- geo[,-1]
sum(rownames(geo) == rownames(bray))


library(ggplot2)
library(ggpubr)#correlation  add p3  相关
library(reshape2)
list <- list()
list2 <- list()
for (i in order) {
  brayd <- bray[rownames(meta)[meta$order == i],rownames(meta)[meta$order == i]]
  phydisd <- phydis[rownames(meta)[meta$order == i],rownames(meta)[meta$order == i]]
  geod <- geo[rownames(meta)[meta$order == i],rownames(meta)[meta$order == i]]
  geocld <- geocl[rownames(meta)[meta$order == i],rownames(meta)[meta$order == i]]
  sum(colnames(phydisd) == colnames(brayd))
  brayd <- as.dist(brayd)
  phydisd <- as.dist(phydisd)
  geod <- as.dist(geod)
  geocld <- as.dist(geocld)
  data <- data.frame(bray = as.vector(brayd),
                     phy = as.vector(phydisd),
                     geo = as.vector(geod),
                     geocl = as.vector(geocld))
  data <- data[data$bray <= 1,]
  data <- na.omit(data)
  if (nrow(data) > 10000){
#    data <- data[sample(1:nrow(data),10000),]
    data <- data[seq(2,nrow(data),floor(nrow(data)/10000)),]
  }
  newdata <- data
  newdata$phy <-  (newdata$phy-min(newdata$phy)) / (max(newdata$phy)-min(newdata$phy))
  newdata$geo <-  (newdata$geo-min(newdata$geo)) / (max(newdata$geo)-min(newdata$geo))
  newdata$geocl <-  (newdata$geocl-min(newdata$geocl)) / (max(newdata$geocl)-min(newdata$geocl))
  colnames(newdata) <- c("Weighted_unifrac_distance","Phylogenetic distance",
                      "Geographical distance","Climatic distance")
  newdata <- melt(newdata, id = c("Weighted_unifrac_distance"))
  p1 <- ggplot(newdata,aes(value,Weighted_unifrac_distance,color = variable))+
    geom_point(size = 1,alpha = 0.1) +# facet_wrap(~variable) +
    labs( x = "Normalized distance", y = "Weighted unifrac distance",
          color = "") +
    geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE) +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'spearman', 
             label.x.npc = 'left', label.y.npc = 'top', size = 5)
  p1
  p <- cor.test(data$bray,data$phy,method = "spearman")
  print(p$estimate)
  list[[i]] <- p1
  data$order <- i
  list2[[i]] <- data
}

library(patchwork)
p <- list[[1]] + list[[2]] + list[[3]] + list[[4]] + 
  list[[5]] + list[[6]] + list[[7]]
p
ggsave(p, file="order2.pdf", width=20, height=12)

#start mixed effect model
data <- list2[[1]]
for (i in 2:7) {
  data <- rbind(data,list2[[i]])
}
data$geo <- scale(data$geo)
colnames(data) <- c("bray","phy","geo","geocl","order")
# data2 <- data.frame(matrix(0,nrow = nrow(data),ncol = ncol(data)))
# colnames(data2) <- c("bray","phy","geo","geocl","order")
# data2[1:nrow(data2),] <- data[1:nrow(data),]
library(lme4)
library(tidyverse)
library(afex)

library(glmm.hp)
model1<-lmer(bray ~  1 + phy +  geocl + geo +
            (1 + phy + geocl + geo |order),data=data)
model3<-lmer(bray ~  1  +  geocl + geo +
               (1 + phy + geocl + geo |order),data=data)
like_test <- anova(model3,model1)

library(effectsize)
result <- effectsize(model1) # need long time
coef(model1)
#result3 <- r.squaredGLMM(rt_full.mod4)
#result2 <- glmm.hp(rt_full.mod4)
summary(model1)
#单因素
#anova(rt_full.mod4,type="I")
result
result <- result[!(result$Parameter == "(Intercept)"),]
result$name <- c("Phylogenetic distance","Climatic distance","Geographical distance")
mixed_p <- ggplot(result, aes(x = name, y = Std_Coefficient, color = name)) +
  geom_point(size = 6) + 
  scale_color_manual(values = c("Phylogenetic distance" = "#F8766D",
                                "Climatic distance" = "#619CFF",
                                "Geographical distance" = "#00BA38")) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.3, size =1 ) +  #误差线表示回归系数的标准误
#  scale_color_manual(values = c("#F39B7FFF","#91D1C2FF"), limits = c('soil', 'micro')) +
  coord_flip() +  #横纵轴转置
  geom_hline(yintercept = 0, linetype = 2,size=1) + 
  labs(x = '', y = 'Effect size', color = '') +  
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line.x = element_line(), axis.ticks.y = element_blank(), 
        legend.key = element_blank(),axis.text.y = element_blank(),
        axis.title=element_text(size=14),axis.text=element_text(size=12),
        legend.text = element_text(size=12)) #+
#  scale_x_discrete(breaks = result$Parameter, labels = result$Parameter, position = 'top') +
#  scale_y_continuous(expand = c(0, 0), 
                     #limits = c(-0.1, 0.8)
#  )#
mixed_p
ggsave(mixed_p, file="mixed.pdf", width=4, height=4)
save(mixed_p,file = "mixed_p.Rdata")#

#
# result4 <-  mixed(bray ~ 1 + phy +  geocl + geo +
#         (1 + phy + geocl + geo |order), 
#       data = data, 
# #      control = lmerControl(optimizer = "bobyqa"), 
#       method = 'LRT')
# #phy < 0.001

model2 <- mixed(bray ~  1 + phy +  geocl + geo +
                  (1 + phy + geocl + geo |order),data=data)

rm(bray,geo,geocl,phydis,brayd,geocld,geod,phydisd)
save.image("mixed_p.Rdata")
