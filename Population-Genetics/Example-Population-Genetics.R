######## Population genetics example ########

#### Required packages

library(MASS)
library(abcrf)
library(ranger)
library(data.table)
ncores <- 40

nReplicate <- 10
resBagging <- resRF <- resLVIRF <- resCsrf50 <-  resCsrf150 <- resCsrf250 <- resDVSRF1 <- resDVSRF2 <- resKVRF1 <- resKVRF2 <- resKVRF3 <- resKVRF4 <- resNNRF1 <- resNNRF2 <- resNNRF3 <- rep(0,nReplicate)

sink("Example-Population-Genetics-Output.txt")

data.all <- fread("data.csv")
data.all$mod <- as.factor(data.all$mod)
n.tot <- dim(data.all)[1]
n <- 10000 
n.remain <- n.tot - n

for (k in 1:nReplicate)  {
  
cat(k)
cat("\n")

set.seed(1974+k)
  
### Training data choice

indicesTrain <- sample(1:n.tot, n, replace=FALSE)
indicesTestRemaining <- c(1:n.tot)[-indicesTrain]

data.train <- data.all[indicesTrain,]
data.test <- data.all[indicesTestRemaining,]

model.lda <- lda(mod~. , data.train)
xtest.lda <- predict(model.lda, data.test)$x

nTest <- 5000
nTestScen1 <- ceiling(nTest/3)
nTestScen2 <- ceiling(nTest/3)
nTestScen3 <- nTest - 2*ceiling(nTest/3)
idxWindows <- c(1:n.remain)[-1<xtest.lda[,1] & xtest.lda[,1]<1 & -1<xtest.lda[,2] & xtest.lda[,2]<1]

indiceTestScen1 <- sample(idxWindows[data.test$mod[idxWindows] == 1], nTestScen1, replace = FALSE)
indiceTestScen2 <- sample(idxWindows[data.test$mod[idxWindows] == 2], nTestScen2, replace = FALSE)
indiceTestScen3 <- sample(idxWindows[data.test$mod[idxWindows] == 3], nTestScen3, replace = FALSE)
indiceTest <- c(indiceTestScen1, indiceTestScen2, indiceTestScen3)

data.test <- data.test[indiceTest,]

x.train <- as.matrix(data.train[,-1])
x.test <- as.matrix(data.test[,-1])
classeTest <- data.test$mod


### Bagging

baggedRf <- ranger(formula = mod~., data = data.train, num.trees = 100, 
                   mtry = dim(x.train)[2], num.threads = ncores)
predBagging <- predict(object = baggedRf, data = data.test, num.threads = ncores)
resBagging[k] <- mean(predBagging$predictions != classeTest)
cat("Bagging ok\n")

#### Random Forests

classicRF <- ranger(formula = mod~., data = data.train,
                    num.trees = 100, num.threads = ncores)
predRF <- predict(object = classicRF, data = data.test, num.threads = ncores)
resRF[k] <- mean(predRF$predictions != classeTest)
cat("Classic RF ok\n")

#### Local variable importance RF

source("LocalVarImpRF.R")

rf.ranger <- ranger(mod ~ ., data = data.train, num.trees = 100, num.threads = ncores)
impxStd <- matrix(NA, nrow = nTest, ncol=dim(x.train)[2])
for (i in 1:nTest) 
  {
  impxStd[i,] <- LocalVarImp(rf.ranger, x.test[i,,drop=FALSE])
  }
predLVIRF <- factor(c(),levels=levels(data.train$mod))
for (i in 1:nTest)
  {
  rf.local.ranger <- ranger(mod ~ ., data = data.train, num.trees = 100,  split.select.weights = impxStd[i,], num.threads = ncores)
  predLVIRF[i] <- predict(rf.local.ranger, data=x.test[i,,drop=FALSE])$predictions
  }
resLVIRF[k] <- mean(predLVIRF != classeTest)
cat("Local variable importance RF ok\n")

#### Case Specific Random Forests

## Nmin = 50

predCsrf50 <- csrf(mod~., training_data = data.train, test_data = data.test, params1 = list(num.trees=100, mtry = 48, min.node.size = 50, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf50[k] <- mean(predCsrf50 != classeTest)

## Nmin = 150

predCsrf150 <- csrf(mod~., training_data = data.train, test_data = data.test, params1 = list(num.trees=100, mtry = 48, min.node.size = 150, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf150[k] <- mean(predCsrf150 != classeTest)

## Nmin = 250

predCsrf250 <- csrf(mod~., training_data = data.train, test_data = data.test, params1 = list(num.trees=100, mtry = 48, min.node.size = 250, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf250[k] <- mean(predCsrf250 != classeTest)

cat("Case Specific Random Forests ok\n")

#### Local dynamic selection RF

source("DynamicVotingWithSelectionRF.R")

## 3000 neighbors, we keep 100 best trees (all)

predDVSRF1 <- dynamicVoting(formula = mod~., data = data.train, dataTest = x.test, K = 3000, ntree = 100, ntreeToKeep = 100, ncores = ncores)
resDVSRF1[k] <- mean(predDVSRF1$prediction !=  classeTest)

## 3000 neighbors, we keep 50 best trees

predDVSRF2 <- dynamicVoting(formula = mod~., data = data.train, dataTest = x.test, K = 3000, ntree = 100, ntreeToKeep = 50, ncores = ncores)
resDVSRF2[k] <- mean(predDVSRF2$prediction !=  classeTest)
cat("Local dynamic selection RF ok\n")

#### Kernel voting

source("KernelVotingRF.R")

## alpha = 1

predKVRF1 <- kernelVoting(formula = mod~., data = data.train, dataTest = x.test, ntree = 100, ncores = ncores, rule = "quantile", alpha = 1)
resKVRF1[k] <- mean(predKVRF1$prediction != classeTest)

## alpha = 0.75

predKVRF2 <- kernelVoting(formula = mod~., data = data.train, dataTest = x.test, ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.75)
resKVRF2[k] <- mean(predKVRF2$prediction != classeTest)

## alpha = 0.5

predKVRF3 <- kernelVoting(formula = mod~., data = data.train, dataTest = x.test, ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.5)
resKVRF3[k] <- mean(predKVRF3$prediction != classeTest)

## alpha = 0.25

predKVRF4 <- kernelVoting(formula = mod~., data = data.train, dataTest = x.test, ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.25)
resKVRF4[k] <- mean(predKVRF4$prediction != classeTest)
cat("Kernel voting ok\n")

#### Nearest-neighbors followed by classic RF

madInit <- apply(X = data.train[,-1], 2, mad)

## 1000 NN

K <- 1000
predNNRF1 <- factor(c(),levels=levels(data.train$mod))
# tp <- txtProgressBar(min = 1, max = nTest, style = 3, char = "*")
for(i in 1:nTest)
  {
  distances <- sapply(1:n, function(X) sqrt(mean(((x.train[X,]-x.test[i,])/madInit)^2)))
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.train[toKeep,]
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  predNNRF1[i] <- predict(rfNN, data=data.frame(x.test[i,,drop=FALSE]), num.threads=ncores)$predictions
#  setTxtProgressBar(tp, i)
  }
resNNRF1[k] <- mean(predNNRF1 != classeTest)

## 1500 NN

K <- 1500
predNNRF2 <- factor(c(),levels=levels(data.train$mod))
# tp <- txtProgressBar(min = 1, max = nTest, style = 3, char = "*")
for(i in 1:nTest)
{
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.train[X,]-x.test[i,])/madInit )^2)) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.train[toKeep,]
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  predNNRF2[i] <- predict(rfNN, data=data.frame(x.test[i,,drop=FALSE]), num.threads=ncores)$predictions
#  setTxtProgressBar(tp, i)
}
resNNRF2[k] <- mean(predNNRF2 != classeTest)

# 2500 NN

K <- 2500
predNNRF3 <- factor(c(),levels=levels(data.train$mod))
# tp <- txtProgressBar(min = 1, max = nTest, style = 3, char = "*")
for(i in 1:nTest)
{
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.train[X,]-x.test[i,])/madInit )^2)) )
  ord <- order(distances)
  toKeep <- ord[1:K]
  data.trainNN <- data.train[toKeep,]
  rfNN <- ranger(formula = mod~., data = data.trainNN, num.trees = 100, num.threads=ncores)
  predNNRF3[i] <- predict(rfNN, data=data.frame(x.test[i,,drop=FALSE]), num.threads=ncores)$predictions
#  setTxtProgressBar(tp, i)
}
resNNRF3[k] <- mean(predNNRF3 != classeTest)
cat("Nearest-neighbors followed by classic RF\n")

write.table(cbind(resBagging,resRF,resLVIRF,resCsrf50,resCsrf150,resCsrf250,resDVSRF1,resDVSRF2,resKVRF1,resKVRF2,resKVRF3,resKVRF4,resNNRF1,resNNRF2,resNNRF3),file="Example-Population-Genetics-Res.txt",quote=FALSE,row.names=FALSE)

}

sink()
