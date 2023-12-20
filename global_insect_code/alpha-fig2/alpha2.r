# Calculate the alpha diversity of the population (100 samples) 
library(vegan)
library(picante)
library(Matrix)
base = exp(1)
alpha <- function(x) {
  x <- otu[x,,drop=F]
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[3, ]
  ACE <- est[5, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')	#Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}
load("ASV_tabs1_own.Rdata")
otuz <- ASV_tabs1_own[[1]]   #load asv table
tree <- read.tree('../../new_15928samples/tree/otu.tree')
#read tree file

# Functions for computing alpha diversity
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[3, ]
  ACE <- est[5, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')	#Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}
otu <- as(otuz, 'matrix')
otu <- as.data.frame(otu)
rowSums(otu)
meta <- read.delim("../../meta.txt",
                   header = T)  #read meta information of samples
meta <- meta[,c("X","order","protein_taxon_id")]
rownames(meta) <- meta$X
meta <- meta[meta$order != "Primates",]
otu <- otu[rownames(meta),colSums(otu) > 0]
#meta <- meta[rownames(otu),]
sum(rownames(meta) == rownames(otu))

order <- unique(meta$order)

dis1 <- read.delim("../../new_15928samples/tree/ML.txt",row.names = 1)
#Read the phylogenetic distance file
colnames(dis1) <- rownames(dis1)
result_whole <- data.frame(matrix(-100,nrow = 100,ncol = 9))
colnames(result_whole)[9] <- "dis" 
for (i in 1:100) {
  set.seed(i)
  n <- sample(1:nrow(meta),100, replace = T)
  otuf <- otu[n,]
  otuf <- otuf[,colSums(otuf) > 0]
  z <- meta[n,]
  sum_dis <- 0
  for (j in 1:99) {
    for (k in (j+1):100) {
      sum_dis1 <- dis1[z$protein_taxon_id[j],z$protein_taxon_id[k]]
      sum_dis <- sum_dis + sum_dis1
    }
  }
#  otuf <- cbind(otuf,z)
#  otuf <- aggregate(.~order, data=otuf ,sum) 
  otuf <- colSums(otuf)
  otuf <- as.matrix(otuf)
  otuf <- t(otuf)
  result_alpha <- alpha(otuf,tree)
  result_whole[i,] <- c(result_alpha,sum_dis)
}
colnames(result_whole) <- c(colnames(result_alpha),"dis")

cor.test(result_whole$Richness,result_whole$dis,method = "spearman")
library(ggplot2)
ggplot(result_whole, aes(x = dis, y =  Richness)) + 
  geom_point()

#order seperate useless
write.table(result_whole,"alpha100.txt",col.names=NA, sep = "\t", quote = F)


meta2 <- meta
otu2 <- otu
for (o in order) {
  meta <- meta2[meta2$order == o,]
  otu <- otu2[meta2$order == o,]
  result_whole <- data.frame(matrix(-100,nrow = 100,ncol = 9))
  colnames(result_whole)[9] <- "dis" 
  for (i in 1:100) {
    set.seed(i)
    n <- sample(1:nrow(meta),100, replace = T)
    otuf <- otu[n,]
    otuf <- otuf[,colSums(otuf) > 0]
    z <- meta[n,]
    sum_dis <- 0
    for (j in 1:99) {
      for (k in (j+1):100) {
        sum_dis1 <- dis1[z$protein_taxon_id[j],z$protein_taxon_id[k]]
        sum_dis <- sum_dis + sum_dis1
      }
    }
    #  otuf <- cbind(otuf,z)
    #  otuf <- aggregate(.~order, data=otuf ,sum) 
    otuf <- colSums(otuf)
    otuf <- as.matrix(otuf)
    otuf <- t(otuf)
    result_alpha <- alpha(otuf,tree)
    result_whole[i,] <- c(result_alpha,sum_dis)
  }
  colnames(result_whole) <- c(colnames(result_alpha),"dis")
  list1[[o]] <- result_whole
  cor.test(result_whole$Richness,result_whole$dis,method = "spearman")
  library(ggplot2)
  ggplot(result_whole, aes(x = dis, y =  Richness)) + 
    geom_point()
}
result_order <- list1[[1]]
write.table(result_whole,"alpha100.txt",col.names=NA, sep = "\t", quote = F)