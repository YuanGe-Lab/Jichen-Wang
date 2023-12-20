# code for cubist model
load("RF.Rdata")
# model construction
myfun1 <- function(x){
#  library(caret)
  data <- meta[names(asv_core[[names(asv_core)[x]]]),-1]
  data$asv <-  asv_core[[names(asv_core)[x]]]
  # set.seed(1)
  # model <- train(
  #   asv ~ .,
  #   data = training,
  #   method = 'cubist'
  # )
  # summary(model)
  ctrl <- trainControl(
    method = "cv",
    number = 10,
  )
  R1 <- 0
  err <- c()
  #1:10
  for (i in 1:10) {
    set.seed(i)
    inTraining <- createDataPartition(data$asv, p = .80, list = FALSE)
    training <- data[inTraining,]
    testing  <- data[-inTraining,]
    tuneGrid <- expand.grid(
      committees = c(5, 100),
      neighbors = c(5, 9) #'neighbors' must be less than 10
    )
    set.seed(1)
    model <- train(
      asv ~ .,
      data = training,
      method = 'cubist',
      tuneGrid = tuneGrid,
      #    preProcess = c("center", "scale"),
      trControl = ctrl
    )
#    summary(model)
    predictions1 = predict(model, newdata = training)
    target1 <- training$asv
    R3 <- R2(target1, predictions1)
    predictions2 = predict(model, newdata = testing)
    target2 <- testing$asv
    R4 <- R2(target2, predictions2)
    R5 <- R3 + R4
    # print(R3)
    # print(R4)
    # print("-----------")
    err <- c(err,R3,R4)
    if (R5 > R1){
      R1 <- R5
      finalmodel <- model
      ftraining <- training
      ftesting <- testing
    }
  }
  predictions1 = predict(finalmodel, newdata = ftraining)
  target1 <- ftraining$asv
  trainR2 <- R2(target1, predictions1)
  predictions2 = predict(finalmodel, newdata = ftesting)
  target2 <- ftesting$asv
  testR2 <- R2(target2, predictions2)
  # RMSE
  trainRMSE <- sqrt(mean((target1 - predictions1)^2))
  testRMSE <- sqrt(mean((target2 - predictions2)^2))
  bestTune <- finalmodel[["bestTune"]]
  result <- list(model = finalmodel,
                 training  = ftraining,
                 testing = ftesting,
                 trainRMSE = trainRMSE,
                 testRMSE = testRMSE,
                 trainR2 = trainR2,
                 testR2 = testR2,
                 err = err)
  # time.start <- Sys.time()
  # predictions = predict(model, newdata = testing2)
  # time.end <- Sys.time()
  # time.running <- time.end - time.start
  # print(time.running)
}

library(snowfall)
sfInit(parallel = T, cpus = 4)
sfExport("meta","asv_core")
sfLibrary(caret)
result <- sfLapply(x = 1:length(asv_core),fun = myfun1)
sfStop()
for (i in 1:length(result)) {
  print(result[[i]]$testR2)
}
###prediction##
world_point <- read.delim("../point.txt",header = T, row.names = 1)
world_point <- world_point[,c("lon","lat",colnames(meta)[2:ncol(meta)])]
rownames(world_point) <- paste0("S",rownames(world_point))
world_point <- na.omit(world_point)
site <- world_point[,c("lon","lat")]
world_point <- world_point[,-c(1:2)]


myfun2 <- function(x){
  library(caret)
  predictions = predict(result[[x]]$model, newdata = world_point)
}
library(snowfall)
sfInit(parallel = T, cpus = 4)
sfExport("result","world_point")
sfLibrary(caret)
result2 <- sfLapply(x = 1:4,fun = myfun2)
sfStop()

save.image("cubist.Rdata")#保存所有
