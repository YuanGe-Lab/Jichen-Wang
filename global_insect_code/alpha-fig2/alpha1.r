# Code for calculating the correlation between alpha diversity distance and genetic diversity distance
alpha <- read.delim("../alpha.txt",row.names = 1)
dis1 <- read.delim("../../new_15928samples/tree/ML.txt",row.names = 1)
colnames(dis1) <- rownames(dis1)
meta <- read.delim("../../meta.txt",header = T)
meta <- meta[,c("X","order","protein_taxon_id")]
rownames(meta) <- meta$X
meta <- meta[meta$order != "Primates",]
alpha <- alpha[rownames(meta),]

result <- data.frame(matrix(-100,nrow = 200,ncol = 10))
colnames(result) <- c("dis",colnames(alpha))
for (i in 1:200) {
  set.seed(i)
  sample1 <- sample(1:nrow(meta),2)
#  print(meta$X[sample1[1]])
#  print(meta$X[sample1[2]])
  dis <- dis1[meta$protein_taxon_id[sample1[1]],meta$protein_taxon_id[sample1[2]]]
  result[i,1] <- dis 
  for (j in 1:ncol(alpha)) {
    if (is.na(alpha[sample1[1],j]) | is.na(alpha[sample1[2],j])){
      euclidean <- -100
    }else{
      euclidean <- abs(alpha[sample1[1],j] - alpha[sample1[2],j])
      result[i,j + 1] <- euclidean
    }
  }
}
myfun <- function(x) {
  return(min(x))
}
result <-  result[apply(result, 1, myfun) != -100,]
library(reshape2)
result <- subset(result, select = -c(goods_coverage,qiime_pd)) 
colnames(result)[colnames(result) == "PD_whole_tree"] <- "Phylogenetic diversity"
result_backups <-result
for (i in 2:ncol(result)) {
  result[,i] <- (result[,i] - min(result[,i]))/(max(result[,i] - min(result[,i])))
}

result <- melt(result, id = c("dis"))

library(ggplot2)
library(ggpubr)#correlation  add p3  相关
result <- result[result$variable %in% c("Richness","ACE","Simpson"),] #only three p < 0.05
p <- ggplot(result, aes(x = dis, y =  value, color = variable)) + 
  geom_point(alpha = 0.5) +
  scale_color_manual(values=c("Richness" = "#66C2A5",
                              "Simpson" = "#FC8D62",
                              "ACE" = "#E78AC3")) + 
  geom_smooth(method = 'lm', formula = y~x, se = F, show.legend = FALSE) +
  labs(x = "Phylogenetic distance", y = "Normalized alpha distance",color = "") +
#  facet_wrap(~variable,scales = "free") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = '~`,`~')), method = 'spearman', 
           label.x.npc = 'left', label.y.npc = 'top', size = 5,show.legend=FALSE,
           r.accuracy = 0.01, p.accuracy = 0.01) +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text=element_text(size=12))
p
ggsave(p, file="alpha_dis.pdf", width=6, height=4)
