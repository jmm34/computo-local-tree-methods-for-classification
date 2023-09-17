######## Unbalanced Gaussian example without Noise ########

#### Required packages

library(mvtnorm)
library(ranger,lib="~/scratch/soft")

args <- commandArgs(trailingOnly = TRUE)
ncores=as.numeric(args[1])


nReplicate <- 10
resBagging <- resRF <- resLVIRF <- resCsrf5 <- resCsrf10 <- resCsrf50 <-  resCsrf150 <- resCsrf250 <- resCsrf350 <- resDVSRF1 <- resDVSRF2 <- resKVRF1 <- resKVRF2 <- resKVRF3 <- resKVRF4 <- resNNRF1 <- resNNRF2 <- resNNRF3 <- rep(0,nReplicate)

timeBagging <- timeRF <- timeLVIRF <- timeCsrf5 <- timeCsrf10 <- timeCsrf50 <-  timeCsrf150 <- timeCsrf250 <- timeCsrf350 <- timeDVSRF1 <- timeDVSRF2 <- timeKVRF1 <- timeKVRF2 <- timeKVRF3 <- timeKVRF4 <- timeNNRF <- rep(0,nReplicate)


for (k in 7:nReplicate)  {
  
cat(k)
cat("\n")

set.seed(1974+k)
  
### Training data generation

num.train <- 10000
n=num.train
x.train <- rmvnorm(num.train,rep(0,3),diag(rep(4,3)))
dist.origin.train <- sqrt(apply(x.train^2,1,sum))
cosx=x.train[,1] / sqrt(apply(x.train[,c(1,3)]^2,1,sum))
cosy=x.train[,2] / sqrt(apply(x.train[,c(1,2)]^2,1,sum))
cosz=x.train[,3] / sqrt(apply(x.train[,c(2,3)]^2,1,sum))
ind3 <- dist.origin.train>3.75
sum(ind3)
ind2 <- dist.origin.train>2.5 & dist.origin.train<=3.75
sum(ind2)
ind1 <- dist.origin.train<=2.5
sum(ind1)

classe.train <- rep(0,num.train)
for (i in 1:num.train)
{
  if (ind1[i] )
  { 
  	if (x.train[i,1]>0 & cosy[i]>(-1/2) ) 
  		{classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.8,0.1,0.1)) } else   	if (x.train[i,1]<0 & cosy[i]>(-1/2) ) 
  		{classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.8,0.1)) } else  classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.1,0.8))
  }	
  if (ind2[i]) 
  { 
  	if (x.train[i,2]>0 & cosz[i]>(-1/2) ) 
  		{classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.8,0.1,0.1)) } else   	if (x.train[i,2]<0 & cosz[i]>(-1/2) ) 
  		{classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.8,0.1)) } else  classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.1,0.8))
  }	 
  if (ind3[i]) 
    { 
  	if (x.train[i,3]>0 & cosx[i]>(-1/2) ) 
  		{classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.8,0.1,0.1)) } else   	if (x.train[i,3]<0 & cosx[i]>(-1/2) ) 
  		{classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.8,0.1)) } else  classe.train[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.1,0.8))
  }	
}
#library(rgl)
#plot3d(x.train[,1],x.train[,2],x.train[,3],col=classe.train)
plot(x.train[,1],x.train[,2],col=classe.train,pch="*")
plot(x.train[,2],x.train[,3],col=classe.train,pch="*")

sum(classe.train==1)
sum(classe.train==2)
sum(classe.train==3)


num.vnoise <- 20
x.trainNoised <- cbind(x.train, matrix(runif(num.train*num.vnoise, 0, 1), nrow=num.train))

num.var <- num.vnoise + 3

x.trainNoised.norma <- x.trainNoised
for (ll in 1:num.var)
{
  norma <- mad(x.trainNoised[,ll])
  x.trainNoised.norma[,ll] <- x.trainNoised[,ll]/norma
  #   x.testNoised[,k] <- x.testNoised[,k]/norma
}

# Define the training data set
training <- data.frame(classe=as.factor(classe.train), x.trainNoised.norma)



### Simulate the test data

num.test <- 5000
nTest <-num.test
x.test <- rmvnorm(num.test,rep(0,3),diag(rep(4,3)))
dist.origin.test <- sqrt(apply(x.test^2,1,sum))
cosx=x.test[,1] / sqrt(apply(x.test[,c(1,3)]^2,1,sum))
cosy=x.test[,2] / sqrt(apply(x.test[,c(1,2)]^2,1,sum))
cosz=x.test[,3] / sqrt(apply(x.test[,c(2,3)]^2,1,sum))
ind3 <- dist.origin.test>3.75
sum(ind3)
ind2 <- dist.origin.test>2.5 & dist.origin.test<=3.75
sum(ind2)
ind1 <- dist.origin.test<=2.5
sum(ind1)

classe.test <- rep(0,num.test)
for (i in 1:num.test)
{
  if (ind1[i] )
  { 
  	if (x.test[i,1]>0 & cosy[i]>(-1/2) ) 
  		{classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.8,0.1,0.1)) } else   	if (x.test[i,1]<0 & cosy[i]>(-1/2) ) 
  		{classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.8,0.1)) } else  classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.1,0.8))
  }	
  if (ind2[i]) 
  { 
  	if (x.test[i,2]>0 & cosz[i]>(-1/2) ) 
  		{classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.8,0.1,0.1)) } else   	if (x.test[i,2]<0 & cosz[i]>(-1/2) ) 
  		{classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.8,0.1)) } else  classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.1,0.8))
  }	 
  if (ind3[i]) 
    { 
  	if (x.test[i,3]>0 & cosx[i]>(-1/2) ) 
  		{classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.8,0.1,0.1)) } else   	if (x.test[i,3]<0 & cosx[i]>(-1/2) ) 
  		{classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.8,0.1)) } else  classe.test[i] <- sample(c(1,2,3), 1, prob=c(0.1,0.1,0.8))
  }	
}
# set.seed(321)
x.testNoised <- cbind(x.test, matrix(runif(num.test*num.vnoise, 0, 1),
                                     nrow=num.test))
x.testNoised.norma <- x.testNoised
for (ll in 1:num.var)
{
  norma <- mad(x.trainNoised[,ll])
  # norma <- mad(x.testNoised[,k])
  x.testNoised.norma[,ll] <- x.testNoised[,ll]/norma
}

test <- data.frame(classe=as.factor(classe.test), x.testNoised.norma)
classeTest=as.factor(classe.test)

############################################################

### Bagging

t = system.time({
baggedRf <- ranger(formula = classe~., data = training, num.trees = 100, 
                   mtry = dim(training)[2]-1, num.threads = ncores)
predBagging <- predict(object = baggedRf, data = test[,-1], num.threads = ncores)
resBagging[k] <- mean(predBagging$predictions != classeTest)
})
timeBagging[k]<-t[3]
cat("Bagging ok\n")

#### Random Forests
t = system.time({
classicRF <- ranger(formula = classe~., data = training,
                    num.trees = 100, num.threads = ncores)
predRF <- predict(object = classicRF, data = test[,-1], num.threads = ncores)
resRF[k] <- mean(predRF$predictions != classeTest)
})
timeRF[k]<-t[3]
cat("Classic RF ok\n")

#### Local variable importance RF

source("/home/cleynena/LocalRF/LocalVarImpRF.R")
t = system.time({
rf.ranger <- ranger(classe ~ ., data = training, num.trees = 100, num.threads = ncores)
impxStd <- matrix(NA, nrow = nTest, ncol=dim(training)[2]-1)
for (i in 1:nTest) 
  {
  impxStd[i,] <- LocalVarImp(rf.ranger, test[i,-1,drop=FALSE])
  }
predLVIRF <- factor(c(),levels=levels(training$classe))
# tp <- txtProgressBar(min = 1, max = nTest, style = 3, char = "*")
for (i in 1:nTest)
  {
  rf.local.ranger <- ranger(classe ~ ., data = training, num.trees = 100,  split.select.weights = impxStd[i,], num.threads = ncores)
  predLVIRF[i] <- predict(rf.local.ranger, data=test[i,-1,drop=FALSE])$predictions
#  setTxtProgressBar(tp, i)
  }
resLVIRF[k] <- mean(predLVIRF != classeTest)
})
timeLVIRF[k]<-t[3]
cat("Local variable importance RF ok\n")

#### Case Specific Random Forests

## Nmin = 5
t = system.time({
predCsrf5 <- csrf(classe~., training_data = training, test_data = data.frame(test[,-1]), params1 = list(num.trees=100, mtry = dim(training)[2]-1, min.node.size = 5, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf5[k] <- mean(predCsrf5 != classeTest)
})
timeCsrf5[k]<-t[3]

## Nmin = 10
t = system.time({
predCsrf10 <- csrf(classe~., training_data = training, test_data = data.frame(test[,-1]), params1 = list(num.trees=100, mtry = dim(training)[2]-1, min.node.size = 10, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf10[k] <- mean(predCsrf10 != classeTest)
})
timeCsrf10[k]<-t[3]

## Nmin = 50
t = system.time({
predCsrf50 <- csrf(classe~., training_data = training, test_data = data.frame(test[,-1]), params1 = list(num.trees=100, mtry = dim(training)[2]-1, min.node.size = 50, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf50[k] <- mean(predCsrf50 != classeTest)
})
timeCsrf50[k]<-t[3]

## Nmin = 150
t = system.time({
predCsrf150 <- csrf(classe~., training_data = training, test_data = data.frame(test[,-1]), params1 = list(num.trees=100, mtry = dim(training)[2]-1, min.node.size = 150, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf150[k] <- mean(predCsrf150 != classeTest)
})
timeCsrf150[k]<-t[3]

## Nmin = 250
t = system.time({
predCsrf250 <- csrf(classe~., training_data = training, test_data = data.frame(test[,-1]), params1 = list(num.trees=100, mtry = dim(training)[2]-1, min.node.size = 250, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf250[k] <- mean(predCsrf250 != classeTest)
})
timeCsrf250[k]<-t[3]

# Nmin = 350
t = system.time({
predCsrf350 <- csrf(classe~., training_data = training, test_data = data.frame(test[,-1]), params1 = list(num.trees=100, mtry = dim(training)[2]-1, min.node.size = 350, num.threads = ncores), params2 = list(num.trees=100, num.threads = ncores))
resCsrf350[k] <- mean(predCsrf350 != classeTest)
})
timeCsrf350[k]<-t[3]
cat("Case Specific Random Forests ok\n")

#### Local dynamic selection RF

source("/home/cleynena/LocalRF/DynamicVotingWithSelectionRF.R")

## 3000 neighbors, we keep 100 best trees (all)
t = system.time({
predDVSRF1 <- dynamicVoting(formula = classe~., data = training, dataTest = data.frame(test[,-1]), K = 3000, ntree = 100, ntreeToKeep = 100, ncores = ncores)
resDVSRF1[k] <- mean(predDVSRF1$prediction !=  classeTest)
})
timeDVSRF1[k]<-t[3]

t = system.time({
predDVSRF2 <- dynamicVoting(formula = classe~., data = training, dataTest = data.frame(test[,-1]), K = 3000, ntree = 100, ntreeToKeep = 50, ncores = ncores)
resDVSRF2[k] <- mean(predDVSRF2$prediction !=  classeTest)
})
timeDVSRF2[k]<-t[3]
cat("Local dynamic selection RF ok\n")

#### Kernel voting

source("/home/cleynena/LocalRF/KernelVotingRF.R")

## alpha = 1
t = system.time({
predKVRF1 <- kernelVoting(formula = classe~., data = training, dataTest = data.frame(test[,-1]), ntree = 100, ncores = ncores, rule = "quantile", alpha = 1)
resKVRF1[k] <- mean(predKVRF1$prediction != classeTest)
})
timeKVRF1[k]<-t[3]

## alpha = 0.75
t = system.time({
predKVRF2 <- kernelVoting(formula = classe~., data = training, dataTest = data.frame(test[,-1]), ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.75)
resKVRF2[k] <- mean(predKVRF2$prediction != classeTest)
})
timeKVRF2[k]<-t[3]

## alpha = 0.5
t = system.time({
predKVRF3 <- kernelVoting(formula = classe~., data = training, dataTest = data.frame(test[,-1]), ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.5)
resKVRF3[k] <- mean(predKVRF3$prediction != classeTest)
})
timeKVRF3[k]<-t[3]

## alpha = 0.25
t = system.time({
predKVRF4 <- kernelVoting(formula = classe~., data = training, dataTest = data.frame(test[,-1]), ntree = 100, ncores = ncores, rule = "quantile", alpha = 0.25)
resKVRF4[k] <- mean(predKVRF4$prediction != classeTest)
})
timeKVRF4[k]<-t[3]
cat("Kernel voting ok\n")


#### Nearest-neighbors followed by classic RF

madInit <- apply(X = training[,-1], 2, mad)

t = system.time({
predNNRF1 <- factor(c(),levels=levels(training$classe))
predNNRF2 <- factor(c(),levels=levels(training$classe))
predNNRF3 <- factor(c(),levels=levels(training$classe))
# tp <- txtProgressBar(min = 1, max = nTest, style = 3, char = "*")
for(i in 1:nTest)
  {
  distances <- sapply(1:n, function(X) sqrt(mean( ( (training[X,-1]-test[i,-1])/madInit )^2)) )
  ord <- order(distances)
  toKeep1 <- ord[1:1000]
  trainingNN1 <- training[toKeep1,]
  toKeep2 <- ord[1:1500]
  trainingNN2 <- training[toKeep2,]
  toKeep3 <- ord[1:2500]
  trainingNN3 <- training[toKeep3,]    
  rfNN <- ranger(formula = classe~., data = trainingNN1, num.trees = 100, num.threads=ncores)
  predNNRF1[i] <- predict(rfNN, data=data.frame(test[i,-1,drop=FALSE]), num.threads=ncores)$predictions
  rfNN <- ranger(formula = classe~., data = trainingNN2, num.trees = 100, num.threads=ncores)
  predNNRF2[i] <- predict(rfNN, data=data.frame(test[i,-1,drop=FALSE]), num.threads=ncores)$predictions
  rfNN <- ranger(formula = classe~., data = trainingNN3, num.trees = 100, num.threads=ncores)
  predNNRF3[i] <- predict(rfNN, data=data.frame(test[i,-1,drop=FALSE]), num.threads=ncores)$predictions  
  }
resNNRF1[k] <- mean(predNNRF1 != classeTest)
resNNRF2[k] <- mean(predNNRF2 != classeTest)
resNNRF3[k] <- mean(predNNRF3 != classeTest)
})
timeNNRF[k]<-t[3]

cat("Nearest-neighbors followed by classic RF OK\n")


write.table(cbind(resBagging,resRF,resLVIRF,resCsrf5,resCsrf10,resCsrf50,resCsrf150,resCsrf250,resCsrf350,resDVSRF1,resDVSRF2,resKVRF1,resKVRF2,resKVRF3,resKVRF4,resNNRF1,resNNRF2,resNNRF3), file=paste("Spherical.txt",sep=""),quote=FALSE,row.names=FALSE)

write.table(cbind(timeBagging,timeRF,timeLVIRF,timeCsrf5,timeCsrf10,timeCsrf50,timeCsrf150,timeCsrf250,timeCsrf350,timeDVSRF1,timeDVSRF2,timeKVRF1,timeKVRF2,timeKVRF3,timeKVRF4,timeNNRF), file=paste("Spherical-time.txt",sep=""),quote=FALSE,row.names=FALSE)
}


