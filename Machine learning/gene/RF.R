####加载程序???
library(MASS) #breast and pima indian data
library(randomForest) #random forests
library(caret) #tune hyper-parameters
data <-  read.table(file = "uniprot.txt",header = T,row.names=1,sep = "\t")
map <- read.table(file = "meta.txt",header = T,row.names=1)
data <- data.frame(t(data))
me <- merge(data, map, by = "row.names", all = TRUE)
row.names(me) <- me$Row.names
me <- me[,-which(colnames(me) == "Row.names")]
#add class and reformat
train <- me
train$class <- factor(train$class)
#----------交叉验证--------------------------------
max = 0
num = 0
for (i in 1:10){
  set.seed(i)
  folds <- createFolds(y=train$class,k=10)
  for(j in 1:10){
    fold_test <- train[folds[[j]],]
    fold_train <- train[-folds[[j]],]
    set.seed(j)
    fold_pre <- randomForest(class ~ ., data = fold_train)
    fold_predict <- predict(fold_pre,type = 'response',newdata=fold_test)
    t <- table(fold_predict, fold_test$class)
    fold_accuracy =  (sum(diag(t))/sum(t))
    print(fold_accuracy)
    if (fold_accuracy > max) {
      max <- fold_accuracy
      bestmodel <-  fold_pre
    }
  }
}
######按重要性进行排序######
varImpPlot(bestmodel, scale = TRUE,
           main = "Variable Importance Plot - PSA Score")
im <- data.frame(importance(bestmodel))
im$genus <- rownames(im)
imo <- im[order(im[,1],decreasing = T),]
train1 <- train[,-ncol(train)]
train1 <- data.frame(t(train1))
me <- merge(train1, imo, by = "row.names", all = TRUE)
me <- me[order(me$MeanDecreaseGini,decreasing = T),]
row.names(me) <- me$Row.names
me <- me[,-which(colnames(me) == "Row.names")]
me <- me[,-which(colnames(me) == "MeanDecreaseGini")]
me <- me[,-which(colnames(me) == "genus")]
data <- data.frame(t(me))
me <- merge(data, map, by = "row.names", all = TRUE)
row.names(me) <- me$Row.names
me <- me[,-which(colnames(me) == "Row.names")]
#add class and reformat
train <- me
train$class <- factor(train$class)
###########前6个基因#####
train1 <- train[,c(1:6,64)]
#----------交叉验证--------------------------------
max = 0
num = 0
set.seed(6)
folds <- createFolds(y=train1$class,k=10)
for(j in 1:10){
    fold_test <- train1[folds[[j]],]
    fold_train <- train1[-folds[[j]],]
    set.seed(j)
    fold_pre <- randomForest(class ~ ., data = fold_train)
    fold_predict <- predict(fold_pre,type = 'response',newdata=fold_test)
    t <- table(fold_predict, fold_test$class)
    fold_accuracy =  (sum(diag(t))/sum(t))
    print(fold_accuracy)
    if (fold_accuracy > max) {
      max <- fold_accuracy
      bestmodel <-  fold_pre
      trainf <- train1[-folds[[j]],]
      testf <- train1[folds[[j]],]
    }
}
set.seed(4)
fold_predict <- predict(bestmodel,type = 'response',newdata=testf)
t <- table(fold_predict, testf$class)
test_accuracy =  (sum(diag(t))/sum(t))


fold_predict <- predict(bestmodel,type = 'response',newdata=trainf)
t <- table(fold_predict, trainf$class)
train_accuracy =  (sum(diag(t))/sum(t))

library(ROCR)
data.probs <- predict(bestmodel, newdata = testf, type = "prob") #save probabilities
pred.test <- prediction(data.probs[,2], testf$class)
perf.test <- performance(pred.test, "tpr", "fpr")

data.probs <- predict(bestmodel, newdata = trainf, type = "prob") #save probabilities
pred.train <- prediction(data.probs[,2], trainf$class)
perf.train <- performance(pred.train, "tpr", "fpr")
######################
pdf(file="myplot.pdf")
plot(perf.train, main = "ROC", col = "#FF0033")
plot(perf.test, col = "#00FF66", add = T)
abline(0,1,col="#660000",lty=2)
legend("bottomright",c(paste0("Train: ","Accuracy = ",round(train_accuracy,3)," AUC = ",
                              round(as.numeric(performance(pred.train, "auc")@y.values), 3)),
                       paste0("Test: ","Accuracy = ",round(test_accuracy,3)," AUC = ",
                              round(as.numeric(performance(pred.test, "auc")@y.values), 3))), 
       col = c("#FF0033","#00FF66"),
       cex = 1,
       lty = 1,lwd = 4)# cex大小调整
#text(0.5,0.5,paste0("Train: ","Accuracy = ",round(train_accuracy,3)," AUC = ",
#                    round(as.numeric(performance(pred.train, "auc")@y.values), 3)))
#text(0.5,0.4,paste0("Test: ","Accuracy = ",round(test_accuracy,3)," AUC = ",
#                    round(as.numeric(performance(pred.test, "auc")@y.values), 3)))
dev.off()
############################################
pdf(file="important.pdf",width=6, height=5)
varImpPlot(bestmodel, scale = TRUE,
           main = "Variable Importance Plot - PSA Score")
dev.off()
importance(bestmodel)
