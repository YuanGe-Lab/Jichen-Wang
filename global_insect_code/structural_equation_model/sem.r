#  Code for structural equation modeling
library(lavaan)
library(data.table)
li <- c("bray_curtis_distance_matrix", 
        "jaccard_distance_matrix",
        "unweighted_unifrac_distance_matrix",
        "weighted_unifrac_distance_matrix")
#with ML bray 0.62, jac 0.72, unweight 0.62, weight 0.44 
#with ge bray 0.05, jac 0.08 unweight 0.03, weight 0.02 
#with cl2 bray -0.02, jac -0.04 unweight -0.08, weight -0.02  
#with clw bray -0.02, jac -0.09 unweight -0.08, weight -0.02  
i <- 4
bray <- fread(paste0("../core-metrics-results/",li[i],"/distance-matrix.tsv"),
              header = T)
bray <- as.data.frame(bray)
rownames(bray) <- bray$V1
bray <- bray[,-1]
brayd <- as.dist(bray)
#meta <- read.delim("../new_15928samples/meta.txt")
#meta <- meta[,c(1,16,30,31)]
#meta <- meta[meta$order != "Primates",]
di <- c("ML_wd.txt","NJ_wd.txt")
j <- 1
phydis <- fread(paste0("../new_15928samples/tree/",di[j]),
                header = T)
phydis <- as.data.frame(phydis)
rownames(phydis) <- colnames(phydis)
sum(colnames(phydis) == colnames(bray))
phydisd <- as.dist(phydis)

geo <- fread(paste0("../map_geodistance/distance2_whole/site_dis.txt"),
               header = T)
geo <- as.data.frame(geo)
rownames(geo) <- geo$V1
geo <- geo[,-1]
sum(rownames(geo) == rownames(bray))
geod <- as.dist(geo)

k <- 3
gl <- c("site_dis.txt","dis_2.txt","dis_whole.txt")
geocl <- fread(paste0("../map_geodistance/distance2_whole/",gl[k]),
             header = T)
geocl <- as.data.frame(geocl)
rownames(geocl) <- geocl$V1
geocl <- geocl[,-1]
sum(rownames(geocl) == rownames(bray))
geocld <- as.dist(geocl)
rm(bray,geo,geocl,phydis)
data <- data.frame(bray = as.vector(brayd), 
                   geo = as.vector(geod),
                   phy =  as.vector(phydisd),
                   geocl =as.vector(geocld))
data <- data[!is.na(data$geocl),]
data <- data[data$bray <= 1,]

#data <- as.data.frame(scale(data))
rm(brayd,geod,phydisd,geocld)
gc()
seq1 <- seq(1,nrow(data),nrow(data)/100000)
data1 <- data[seq1,]
data1 <- scale(data1)
cor.test(data1[,1],data1[,2],method = 'pearson')
#print(object.size(data), unit="GB") 
model <- ' bray  ~ phy +  geo + geocl   
           phy ~ geo + geocl
           geocl ~ geo'
path <- sem(model = model, data = data1,se = 'bootstrap', bootstrap = 1000)
#cor.test(data$bray,data$phy, method = 'pearson')
library(semPlot)
pdf(file = "whole_new.pdf",
    width = 7,           # 宽
    height = 7,          # 高
    bg = "white",          # 背景颜色
)
semPaths(path, what = 'path', whatLabels = 'stand',
         layout = 'tree', style = 'lisrel', residuals = F, edge.label.cex = 1)
#p <- semPaths(path, what = 'path', whatLabels = 'stand',layout = 'tree', style = 'lisrel', residuals = T, edge.label.cex = 1)
#library(semptools)
#p2 <- mark_sig(p, path)
#plot(p2)
dev.off()
summary(path, fit.measures = TRUE,standardize = T)
parameterestimates(path)
para1 <- fitMeasures(path, c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))
para1 <- data.frame(para1)

data2 <- data[data$phy < 0.07,]
#data2 <- data[data$phy < 0.8,] #no use
seq2 <- seq(1,nrow(data2),nrow(data2)/100000)
data2 <- data2[seq2,]
data2 <- scale(data2)
model <- ' bray  ~ phy +  geo + geocl   
           phy ~ geo + geocl
           geocl ~ geo'
path <- sem(model = model, data = data2,se = 'bootstrap', bootstrap = 1000)
#cor.test(data$bray,data$phy, method = 'pearson')
#library(semPlot)
pdf(file = "host_new.pdf",
    width = 7,           # 宽
    height = 7,          # 高
    bg = "white",          # 背景颜色
)
semPaths(path, what = 'path', whatLabels = 'stand',
         layout = 'tree', style = 'lisrel', residuals = F, edge.label.cex = 1)
dev.off()
summary(path, fit.measures = TRUE,standardize = T)
#bray  ~   the p value greater than 0.05
parameterestimates(path)
para2 <- fitMeasures(path, c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))
para2 <- data.frame(para2)


#data3 <- data[data$geo < 11500,]
data3 <- data[data$geo < 400000,]
seq3 <- seq(1,nrow(data3),nrow(data3)/100000)
data3 <- data3[seq3,]
data3 <- scale(data3)
model <- ' bray  ~ phy +  geo + geocl   
           phy ~ geo + geocl
           geocl ~ geo'
path <- sem(model = model, data = data3,se = 'bootstrap', bootstrap = 1000)
#cor.test(data$bray,data$phy, method = 'pearson')
library(semPlot)
pdf(file = "geo_new.pdf",
    width = 7,           # 宽
    height = 7,          # 高
    bg = "white",          # 背景颜色
)
semPaths(path, what = 'path', whatLabels = 'stand',
         layout = 'tree', style = 'lisrel', residuals = F, edge.label.cex = 1)
dev.off()
summary(path, fit.measures = TRUE,standardize = T)
# phy  ~ geocl  the p value  0.103
parameterestimates(path)
para3 <- fitMeasures(path, c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))
para3 <- data.frame(para3)

data4 <- data[data$geocl < 0.06,]
seq4 <- seq(1,nrow(data4),nrow(data4)/100000)
data4 <- data4[seq4,]
data4 <- scale(data4)
model <- ' bray  ~ phy +  geo + geocl   
           phy ~ geo + geocl
           geocl ~ geo'
path <- sem(model = model, data = data4,se = 'bootstrap', bootstrap = 1000)
#cor.test(data$bray,data$phy, method = 'pearson')
library(semPlot)
pdf(file = "geocl_new.pdf",
    width = 7,           # 宽
    height = 7,          # 高
    bg = "white",          # 背景颜色
)
semPaths(path, what = 'path', whatLabels = 'stand',
         layout = 'tree', style = 'lisrel', residuals = F, edge.label.cex = 1)
dev.off()
summary(path, fit.measures = TRUE,standardize = T)
#bray  ~   the p value greater than 0.05
parameterestimates(path)
para4 <- fitMeasures(path, c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))
para4 <- data.frame(para4)
para <- cbind(para1,para2,para3,para4)
colnames(para) <- c("Whole","host","geo","geocl")
write.table(para,"parameter.txt",col.names=NA, sep = "\t", quote = F) 
