#define core asv
library(Matrix)
load("../ASV_tabs1_own.Rdata")
tax <- read.delim("classification.txt")
meta <- read.delim("../meta.txt")
unique(meta$class)

#meta2 <- meta[meta$class == "Insecta",]


rownames(ASV_tabs1_own$V4)
rowSums(ASV_tabs1_own$V4)
asv <- as.matrix(ASV_tabs1_own$V4)
ncol(asv)
asv <- asv[,colSums(asv) > 0]
asv <- as.data.frame(asv)
ncol(asv)
#asv_nohuman <- asv[meta2$X,]
order <- unique(meta$order)
order2 <- as.data.frame(table(meta$order))
species <- list()
nu <- 0
for (i in order) {
  meta2 <- meta[meta$order == i,]
  species2 <- as.data.frame(table(meta2$species))
  species2 <- species2[(species2$Freq > 9) & (species2$Var1 != ""),]
  species2$Var1 <- as.character(species2$Var1)
  species[[i]] <- species2
  nu <- nu +  sum(species2$Freq)
}
print(nu)
# nu = 14409 - human 3446, end 10963 (119 sepcies)

list <- list()
for (i in order) {
  list[[i]] <- list()
  for (j in species[[i]]$Var1) {
    result <- asv[meta$X[meta$species == j],]
    result <- result[,colSums(result) > 0]
    result2 <- result
    o <- data.frame(asv = colnames(result2),
                    sum = colSums(result2))
    o <-  o[order(-o$sum),]
    oli <-  o$asv[1:ceiling(nrow(o) / 5)]
    result[result > 0] <- 1
    result3 <- result2[,colSums(result) > (nrow(result) / 3)]
    result3 <- result3[,colnames(result3) %in% oli]
    list[[i]][[j]] <- result3 
  }
}
save(list,file="species_core.Rdata")


list2 <- c()
nu_w <- c()
nu_t <- 0
for (i in order) {
  nu <- c()
  for (j in species[[i]]$Var1) {
    for (k in colnames(list[[i]][[j]])){
      if (k %in% names(nu)){
        nu[k] <- nu[k] + 1
      }else{
        nu[k] <- 1
      }
      if (i != order[8]){
        if (k %in% names(nu_w)){
          nu_w[k] <- nu_w[k] + 1
        }else{
          nu_w[k] <- 1
        }
      }
    }
  }
  nu <- nu[nu >= (length(species[[i]]$Var1) / 3)]
  list2[[i]] <- nu
  nu_t <- nu_t + nrow(species[[i]])
}

nu_w <- nu_w[nu_w >= ((nu_t - 1) / 3)]
# V4.00007 is the core asv for whole insect
#sample number
nu1 <- 0
for (i in species) {
  n <-  sum(i$Freq)
  nu1 <- nu1 + n
  print(nu1)
}
#spceis number
nu1 <- 0
for (i in species) {
  n <-  nrow(i)
  nu1 <- nu1 + n
  print(nu1)
}
number2 <- data.frame(matrix(1,nrow = nu1,ncol = 3))
colnames(number2) <- c("species","core","order")
n1 <- 1
for (i in 1:8) {
  for (j in 1:length(list[[order[i]]])) {
    number2[n1, 1]  <- names(list[[order[i]]])[j]
    number2[n1, 2]  <- ncol(list[[order[i]]][[j]])
    number2[n1, 3]  <- order[i]
    n1 <- n1 + 1
  } 
}
number2$seq <- 1:nrow(number2)
library(ggplot2)
number2$order <- factor(number2$order,levels = order)
p <- ggplot(number2,aes(seq,core,fill = order))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("Blattodea" = "#66C2A5","Coleoptera" = "#FC8D62",
                              "Diptera" = "#8DA0CB", "Hemiptera" = "#E78AC3",
                              "Hymenoptera" = "#A6D854", "Lepidoptera" = "#FFD92F",
                              "Orthoptera" =  "#E5C494", "Primates" = "#B3B3B3")) +
  labs(x = "Species", y = "Core asv", fill = "Order") + 
  xlim(0,120)
p


save(p,file="p.Rdata")
ggsave(p, file="species_core.pdf", width=5, height=4)
number1 <-  data.frame(matrix(1,nrow = 8,ncol = 4))
rownames(number1) <- order
colnames(number1) <- c("species","total","mean","order_core")
for (i in 1:8) {
  n1 <-  length(list[[order[i]]])
  print(n1)
  tax1 <- c()
  n3 <- 0
  for (j in list[[order[i]]]) {
    tax1 <- c(tax1,colnames(j))
    tax1 <- unique(tax1)
    n2 <- length(tax1)
    n3 <- n3 + ncol(j)
  }
  n3 <- (n3 / n1)
  n4 <- length(list2[[order[i]]])
  number1[i,] <- c(n1,n2,n3,n4)
}
write.table(number1,"number.txt",col.names=NA, sep = "\t", quote = F)
coretax <- list2
save(coretax,file="coretax.Rdata")
tax <- c()
for (i in order){
  result4 <-  names(coretax[[i]])
  tax <- c(tax,result4)
}
length(unique(tax))
asv <- asv[,colnames(asv) %in% tax]
asv <- asv[rowSums(asv) > 0,colSums(asv) > 0]
asv <- data.frame(t(asv))
write.table(data.frame(ID=rownames(asv),asv),"whole/asv.txt", row.names=F, quote = F, sep = "\t")

library(treeio)
library(ape)
tree <- read.newick("../otu.tree") 
tree <- keep.tip(tree,rownames(asv)) #保留节点
write.tree(tree,"whole/core.tree")
########no human
human <- meta$X[meta$species == "Homo sapiens"]

tax <- c()
for (i in order[1:7]){
  result4 <-  names(coretax[[i]])
  tax <- c(tax,result4)
}
length(unique(tax))

asv <- asv[rownames(asv) %in% tax,!(colnames(asv) %in% human)]
asv <- asv[rowSums(asv) > 0,colSums(asv) > 0]
write.table(data.frame(ID=rownames(asv),asv),"asv.txt", row.names=F, quote = F, sep = "\t")
tree <- keep.tip(tree,rownames(asv)) #保留节点
write.tree(tree,"core.tree")
