#the code for fig. 2d

library(data.table)
li <- c("weighted_unifrac_distance_matrix")
#with ML bray 0.62, jac 0.72, unweight 0.62, weight 0.44 
#with ge bray 0.05, jac 0.08 unweight 0.03, weight 0.02 
#with cl2 bray -0.02, jac -0.04 unweight -0.08, weight -0.02  
#with clw bray -0.02, jac -0.09 unweight -0.08, weight -0.02  
i <- 1
bray <- fread(paste0("../../core-metrics-results/",li[i],"/distance-matrix.tsv"),
              header = T)
bray <- as.data.frame(bray)
rownames(bray) <- bray$V1
bray <- bray[,-1]
sum(rownames(bray) == colnames(bray))
brayd <- as.dist(bray)
#meta <- read.delim("../new_15928samples/meta.txt")
#meta <- meta[,c(1,16,30,31)]
#meta <- meta[meta$order != "Primates",]
di <- c("ML_wd.txt","NJ_wd.txt")
j <- 1
phydis <- fread(paste0("../../new_15928samples/tree/",di[j]),
                header = T)
phydis <- as.data.frame(phydis)
rownames(phydis) <- colnames(phydis)
sum(colnames(phydis) == colnames(bray))
phydisd <- as.dist(phydis)
#cor.test(brayd,phydisd,method = "spearman")
choice <- seq(1,length(brayd),100)
p <- cor.test(brayd[choice],phydisd[choice],method = "spearman")
p$estimate

#length(brayd)
#p <- cor.test(brayd[1:1000000],phydisd[1:1000000],method = "pearson")
#p$estimate

k <- 3
gl <- c("site_dis.txt","dis_2.txt","dis_whole.txt")
geocl <- fread(paste0("../../map_geodistance/distance2_whole/",gl[k]),
             header = T)
geocl <- as.data.frame(geocl)
rownames(geocl) <- geocl$V1
geocl <- geocl[,-1]
sum(rownames(geocl) == rownames(bray))
geocld <- as.dist(geocl)
p <- cor.test(brayd[choice ],geocld[choice ],method = "spearman")
p$estimate
#p <- cor.test(brayd[choice ],geocld[choice ],method = "pearson")
#p$estimate

k <- 1
geo <- fread(paste0("../../map_geodistance/distance2_whole/",gl[k]),
               header = T)
geo <- as.data.frame(geo)
rownames(geo) <- geo$V1
geo <- geo[,-1]
sum(rownames(geo) == rownames(bray))
geod <- as.dist(geo)
p <- cor.test(brayd[choice ],geod[choice ],method = "spearman")
p$estimate

data <- data[data$bray <= 1,]
data$geo[data$geo == 0] <- 1
data <- na.omit(data)
library(ggplot2)
library(ggpubr)#correlation  add p3  相关
# p1 <- ggplot(data,aes(phy,bray))+
#   geom_point(size = 1) + 
#   labs( x = "Phylogeny distance", y = "Weighted Unifrac") +
#   geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE) +
#   stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'spearman', 
#            label.x.npc = 'left', label.y.npc = 'top', size = 5)
# p1
# 
# p2 <- ggplot(data,aes(log10(geo), bray))+
#   geom_point(size = 1) + 
#   xlim(c(0,NA)) +
#   labs( x = "Log10(geodistance)", y = "Weighted unifrac") +
#   geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE) +
#   stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'spearman', 
#            label.x.npc = 'left', label.y.npc = 'top', size = 5)
# p2
# p3 <- ggplot(data,aes(geocl, bray))+
#   geom_point(size = 1) + 
#     labs( x = "Climate distance", y = "Weighted unifrac") +
#   geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE) +
#   stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'spearman', 
#            label.x.npc = 'left', label.y.npc = 'top', size = 5)
# p3
# library(patchwork)
# p <- p1 + p2 + p3
# p
# ggsave(p, file="whole_new.pdf", width=16, height=4)

data$phy <-  (data$phy-min(data$phy)) / (max(data$phy)-min(data$phy))
data$geo <-  (data$geo-min(data$geo)) / (max(data$geo)-min(data$geo))
data$geocl <-  (data$geocl-min(data$geocl)) / (max(data$geocl)-min(data$geocl))

#data$geo <- scale(data$geo)
#data$geocl <- scale(data$geocl)

library(reshape2)
colnames(data) <- c("Weighted_unifrac_distance","Phylogenetic distance",
                    "Geographical distance","Climatic distance")
data <- melt(data, id = c("Weighted_unifrac_distance"))
p4 <- ggplot(data,aes(value,Weighted_unifrac_distance,color = variable))+
  geom_point(size = 1,alpha = 0.1) +# facet_wrap(~variable) +
  labs( x = "Normalized distance", y = "Weighted unifrac distance",
        color = "") +
  geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'spearman', 
           label.x.npc = 'left', label.y.npc = 'top', size = 5, show.legend=FALSE,
           r.accuracy = 0.01, p.accuracy = 0.01) +
  theme(axis.title=element_text(size=14),axis.text=element_text(size=12),
        legend.text = element_text(size=12))
p4
ggsave(p4, file="whole_new.pdf", width=6, height=4)
# data2 <- data.frame(bray = as.numeric(brayd),
#                    phy =  as.numeric(phydisd),
#                    geo =  as.numeric(geod),
#                    geocl =as.numeric(geocld))
# data2 <- data2[!is.na(rowSums(data2)),]

#library(ppcor)
#partial1 <- pcor.test(x = data2$bray, y = data2$phy, z = data2[,c("geo","geocl")], method = 'spearman')  #partial correlation
#0.348
#partial2 <- pcor.test(x = data2$bray, y = data2$geo, z = data2[,c("phy","geocl")], method = 'spearman')  #partial correlation
#-0.005
#partial3 <- pcor.test(x = data2$bray, y = data2$geocl, z = data2[,c("geo","phy")], method = 'spearman')  #partial correlation
#0.03
