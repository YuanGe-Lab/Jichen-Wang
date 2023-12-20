# code for random_forest
meta <- read.delim("../../../meta.txt",row.names = 1)
meta <- subset(meta,select = -c(dis,Use_intensity,ecosystem2,population))
#meta <- subset(meta,select = -c(dis))
# meta$Use_intensity[meta$Use_intensity == "Light use"] <- 1
# meta$Use_intensity[meta$Use_intensity == "Minimal use"] <- 2
# meta$Use_intensity[meta$Use_intensity == "Intense use"] <- 3
# meta$Use_intensity[meta$Use_intensity == "Cannot decide"] <- NA

meta <- meta[meta$class == "Insecta",]
order <- unique(meta$order)
#meta <- meta[,c(15,31:62)]
meta <- meta[,c(15,31:59)]
meta <- na.omit(meta)
load("../../coretxt.Rdata")
core <- core[,-8]
myfun3 <- function(x){
  sum(x == 0)
}
core <- core[apply(core,1,myfun3) < 7,]

load("../../../ASV_tabs1_own.Rdata")
rownames(ASV_tabs1_own$V4)
rowSums(ASV_tabs1_own$V4)
asv <- as.matrix(ASV_tabs1_own$V4)
ncol(asv)
asv_core <- list()
for (i in 1:nrow(core)) {
  tax <- core[i,][core[i,] != 0]
  sample1 <- rownames(meta)[meta$order %in% tax]
  if (length(sample1) > 500){
    asv_core[[rownames(core)[i]]] <- asv[sample1,rownames(core)[i]]
  }
}
rm(asv,ASV_tabs1_own)
gc()

temp <- list(asv_core = asv_core,
             meta = meta)
save(temp,file="define_environment.Rdata")

myfun1 <- function(){
  library(snowfall)
  sfInit(parallel = T, cpus = 40)
  myfun <- function(x){
    meta2 <- meta[names(asv_core[[names(asv_core)[x]]]),]
    meta2$asv <- asv_core[[names(asv_core)[x]]]
    meta2 <- meta2[,-1]
    set.seed(123)
    forest <- randomForest(asv~., data = meta2, importance = TRUE)
    Variance <- round (100 * forest$rsq[length (forest$rsq)], digits = 2)
    im <- importance(forest)
    if (sum(rownames(im) ==colnames(meta2[,-ncol(meta2)])) == (ncol(meta2) - 1)){
      state <- "right"
    }else{
      state <- "wrong"
    }
    re <- c(names(asv_core)[x],Variance,im[,1],state)
  }
  sfExport("meta","asv_core")
  sfLibrary(randomForest)
  result <- sfLapply(x = 1:length(asv_core),fun = myfun) #列表
  sfStop()
  
  number1 <- data.frame(matrix(0,nrow = length(asv_core),ncol = 2 + ncol(meta)))
  colnames(number1) <- c("asv","variance",colnames(meta)[2:ncol(meta)],"state")
  for (i in 1:length(result)) {
    number1[i,] <- result[[i]]
  }
  sum(number1$state != "right")
  number1 <- number1[number1$variance != "NaN", ]
  number1$variance <- as.numeric(number1$variance)
  number1 <- number1[number1$variance > 30 | number1$asv == "V4.00007",]
  
  im <- data.frame(t(number1[,3:(1+ncol(meta))]))
  colnames(im) <- number1$asv
  for (i in 1:ncol(im)) {
    im[,i] <- as.numeric(im[,i])
  }
  im2 <- im
  for (i in 1:ncol(im)) {
    im <- im[order(im[,i],decreasing = T),]
    im[,i] <- 1:nrow(im)
  }
  
  first <- apply(im, 1,mean)
  second <- apply(im2, 1, mean)
  first <- as.data.frame(first)
  second <- as.data.frame(second)
  first <- first[order(first$first),,drop = F]
  list1 <- list(number1 = number1,
                first = first,
                im = im,
                im2 = im2)
  return(list1)
}

# result2 <- myfun1()
# meta3 <- meta
# #asv <- asv[,result2$number1$asv]
# meta <- meta[,c("order",rownames(result2$first)[1:12])]
load("../define_environment/example.Rdata")
meta <- meta[,c("order",environ)]
result2 <- myfun1()

#normal - no paralell
# library(randomForest)
# for (i in 1:ncol(asv)) {
#   data <- cbind(meta,asv[,i])
#   colnames(data)[ncol(data)] <- "asv"
#   forest <- randomForest(asv~., data = data, importance = TRUE)
#   Variance <- round (100 * forest$rsq[length (forest$rsq)], digits = 2)
#   number1$variance[i] <- Variance
#   im <- importance(forest)
#   number1[i,] <- c(colnames(asv)[i],Variance,im[,1])
#   if (sum(rownames(im) ==colnames(meta)) != ncol(meta)){
#     print("wrong")
#   }
# }
rm(meta3,tax,i,sample1,myfun1,myfun3,temp)
asv_core <- asv_core[result2$number1$asv]
save.image("RF.Rdata")#保存所有
