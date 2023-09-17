######## Unbalanced Gaussian example with Noise ########

#### Required packages

library(mvtnorm)
library(ranger)
library(parallel)
ncores <- detectCores()

nReplicate <- 10
resBayes <- resBagging <- resRF <- resLVIRF <- resCsrf5 <- resCsrf10 <- resCsrf50 <-  resCsrf150 <- resCsrf250 <- resCsrf350 <- resDVSRF1 <- resDVSRF2 <- resKVRF1 <- resKVRF2 <- resKVRF3 <- resKVRF4 <- resNNRF1 <- resNNRF2 <- resNNRF3 <- rep(0,nReplicate)

sink("Example-Gaussian-Unbalanced-With-Noise-Output-10000.txt")

for (k in 1:nReplicate)  {
  
cat(k)
cat("\n")

set.seed(1974+k)
  
### Training data generation

n <- 10000
nNoise <- 100
pi0 <- 0.40
pi1 <- 0.40
pi2 <- 0.10
pi3 <- 0.10
l <- 20
mu0 <- c(c(0.8,3), rep(c(1,2.5), l/2-1))
Sigma0 <- diag(c(c(2,1), rep(c(3,1),l/2-1)) )
mu1 <- c(c(3.2,3), rep(c(2.5,2.5), l/2-1))
Sigma1 <- diag(c(c(2,1), rep(c(3,5),l/2-1)) )
mu2 <- c(c(2,1), rep(c(2,2.3),l/2-1))
Sigma2 <- diag( (rep(c(4,1),l/2) ) )
mu3 <- c(c(2,0), rep(c(2,1.8),l/2-1))
Sigma3 <- diag( (rep(c(2.5,1),l/2) ) )
classe <- sample(x = c(0,1,2,3), size = n, replace = TRUE, prob = c(pi0,pi1,pi2,pi3))
classe <- sort(classe)
n0 <- sum(classe==0)
n1 <- sum(classe==1)
n2 <- sum(classe==2)
n3 <- sum(classe==3)
x.train <- rbind(rmvnorm(n0, mu0, Sigma0), rmvnorm(n1, mu1, Sigma1), rmvnorm(n2, mu2, Sigma2), rmvnorm(n3, mu3, Sigma3))
x.train <- cbind(x.train, matrix(runif(n*nNoise, 0, 1), nrow=n))

nTest <- 5000
classeTest <- sample(c(2,3), size=nTest, prob=c(pi2,pi3), replace=TRUE)
classeTest <- sort(classeTest)

nTest2 <- sum(classeTest==2)
nTest3 <- sum(classeTest==3)

x.test <- rbind(rmvnorm(nTest2, mu2, Sigma2),
                rmvnorm(nTest3, mu3, Sigma3))
x.test <- cbind(x.test, matrix(runif(nTest*nNoise, 0, 1), nrow=nTest))

data.train <- data.frame(mod = as.factor(classe), x.train)
colnames(x.test) <- colnames(data.train)[-1]

#### Bayes classifier

BayesClassifieur <- function(x){
  c0 <- pi0*dmvnorm(x[1:l],mean=mu0,sigma=Sigma0)
  c1 <- pi1*dmvnorm(x[1:l],mean=mu1,sigma=Sigma1)
  c2 <- pi2*dmvnorm(x[1:l],mean=mu2,sigma=Sigma2)
  c3 <- pi3*dmvnorm(x[1:l],mean=mu3,sigma=Sigma3)
  return(c(0,1,2,3)[which.max(c(c0,c1,c2,c3))])
}
predBayes <- rep(NA, nTest)
for(i in 1:nTest) predBayes[i] <- BayesClassifieur(x.test[i,])
resBayes[k] <- mean(predBayes != classeTest)
cat("Bayes ok\n")

### Bagging

baggedRf <- ranger(formula = mod~., data = data.train, num.trees = 100, 
                   mtry = dim(x.train)[2], num.threads = ncores)
predBagging <- predict(object = baggedRf, data = x.test, num.threads = ncores)
resBagging[k] <- mean(predBagging$predictions != classeTest)
cat("Bagging ok\n")

#### Random Forests

classicRF <- ranger(formula = mod~., data = data.train,
                    num.trees = 100, num.threads = ncores)
predRF <- predict(object = classicRF, data = x.test, num.threads = ncores)
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
# tp <- txtProgressBar(min = 1, max = nTest, style = 3, char = "*")
for (i in 1:nTest)
  {
  rf.local.ranger <- ranger(mod ~ ., data = data.train, num.trees = 100,  split.select.weights = impxStd[i,], num.threads = ncores)
  predLVIRF[i] <- predict(rf.local.ranger, data=x.test[i,,drop=FALSE])$predictions
#  setTxtProgressBar(tp, i)
  }
resLVIRF[k] <- mean(predLVIRF != classeTest)
cat("Local variable importance RF ok\n")

#### Case Specific Random Forests

## Nmin = 5

predCsrf5 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test), params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 5, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf5[k] <- mean(predCsrf5 != classeTest)

## Nmin = 10

predCsrf10 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test), params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 10, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf10[k] <- mean(predCsrf10 != classeTest)

## Nmin = 50

predCsrf50 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test), params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 50, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf50[k] <- mean(predCsrf50 != classeTest)

## Nmin = 150

predCsrf150 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test), params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 150, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf150[k] <- mean(predCsrf150 != classeTest)

## Nmin = 250

predCsrf250 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test), params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 250, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf250[k] <- mean(predCsrf250 != classeTest)

# Nmin = 350

predCsrf350 <- csrf(mod~., training_data = data.train, test_data = data.frame(x.test), params1 = list(num.trees=100, mtry = dim(x.train)[2], min.node.size = 350, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf350[k] <- mean(predCsrf350 != classeTest)
cat("Case Specific Random Forests ok\n")

#### Local dynamic selection RF

source("DynamicVotingWithSelectionRF.R")

## 3000 neighbors, we keep 100 best trees (all)

predDVSRF1 <- dynamicVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), K = 3000, ntree = 100, ntreeToKeep = 100, ncores = ncores)
resDVSRF1[k] <- mean(predDVSRF1$prediction !=  classeTest)

## 3000 neighbors, we keep 50 best trees

predDVSRF2 <- dynamicVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), K = 3000, ntree = 100, ntreeToKeep = 50, ncores = ncores)
resDVSRF2[k] <- mean(predDVSRF2$prediction !=  classeTest)
cat("Local dynamic selection RF ok\n")

#### Kernel voting

source("KernelVotingRF.R")

## alpha = 1

predKVRF1 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), ntree = 100, ncores = ncores, rule = "quantile", alpha = 1)
resKVRF1[k] <- mean(predKVRF1$prediction != classeTest)

## alpha = 0.75

predKVRF2 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.75)
resKVRF2[k] <- mean(predKVRF2$prediction != classeTest)

## alpha = 0.5

predKVRF3 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.5)
resKVRF3[k] <- mean(predKVRF3$prediction != classeTest)

## alpha = 0.25

predKVRF4 <- kernelVoting(formula = mod~., data = data.train, dataTest = data.frame(x.test), ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.25)
resKVRF4[k] <- mean(predKVRF4$prediction != classeTest)
cat("Kernel voting ok\n")

#### Nearest-neighbors followed by classic RF

madInit <- apply(X = x.train, 2, mad)

## 1000 NN

K <- 1000
predNNRF1 <- factor(c(),levels=levels(data.train$mod))
# tp <- txtProgressBar(min = 1, max = nTest, style = 3, char = "*")
for(i in 1:nTest)
  {
  distances <- sapply(1:n, function(X) sqrt(mean( ( (x.train[X,]-x.test[i,])/madInit )^2)) )
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

write.table(cbind(resBayes,resBagging,resRF,resLVIRF,resCsrf5,resCsrf10,resCsrf50,resCsrf150,resCsrf250,resCsrf350,resDVSRF1,resDVSRF2,resKVRF1,resKVRF2,resKVRF3,resKVRF4,resNNRF1,resNNRF2,resNNRF3),file="Example-Gaussian-Unbalanced-With-Noise-Res-10000.txt",quote=FALSE,row.names=FALSE)

}

sink()
