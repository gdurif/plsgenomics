### load plsgenomics library
library(plsgenomics)

### generating data
n <- 100
p <- 500
sample1 <- sample.cont(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, 
                       mean.H=0.2, sigma.H=10, sigma.F=5, sigma.E=5)

X <- sample1$X
Y <- sample1$Y

### splitting between learning and testing set
index.train <- sort(sample(1:n, size=round(0.7*n)))
index.test <- (1:n)[-index.train]

Xtrain <- X[index.train,]
Ytrain <- Y[index.train,]

Xtest <- X[index.test,]
Ytest <- Y[index.test,]

### tuning the hyper-parameters
# /!\ on 10 cores
cv1 <- spls.cv(X=Xtrain, Y=Ytrain, lambda.l1.range=seq(0.05, 0.95, by=0.1), ncomp.range=1:5,
               adapt=TRUE, return.grid=TRUE, ncores=10, nfolds=10, nrun=1,
               verbose=FALSE)
str(cv1)

### otpimal values
lambda.l1 <- cv1$lambda.l1.opt
ncomp <- cv1$ncomp.opt

### fitting the model, and predicting new observations
model1 <- spls(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, lambda.l1=lambda.l1, ncomp=ncomp)
str(model1)

### plotting the estimation versus real values for the non centered response
plot(model1$Ytrain, model1$hatY.nc, xlab="real Ytrain", ylab="Ytrain estimates")
points(-1000:1000,-1000:1000, type="l")

### plotting residuals versus centered response values
plot(model1$sYtrain, model1$residuals, xlab="sYtrain", ylab="residuals")

### plotting the predictor coefficients
plot(model1$betahat.nc, xlab="variable index", ylab="coeff")

### mean squares error of prediction on test sample
sYtest <- as.matrix(scale(Ytest, center=model1$meanYtrain, scale=model1$sigmaYtrain))
sum((model1$hatYtest - sYtest)^2) / length(index.test)

### plotting predicted values versus non centered real response values on the test set
plot(model1$hatYtest, sYtest, xlab="real Ytest", ylab="predicted values")
points(-1000:1000,-1000:1000, type="l")

### selected variables
model1$A


######### stability selection procedure
## (/!\ using 8 cores)
stab1 <- spls.stab(X=Xtrain, Y=Ytrain, 
                   lambda.l1.range=seq(0.05, 0.95, by=0.1), ncomp.range=1:5, 
                   ncores=8, nresamp=100, 
                   seed=NULL, verbose=TRUE)
str(stab1)

### heatmap of estimated probabilities
stability.selection.heatmap(stab1)

### selected covariates
tmp <- stability.selection(stab1, piThreshold=0.75, rhoError=10)
tmp

### "true" pertinent covariates vs selected ones
sample1$sel
tmp$selected.predictors

# effect of probability threshold
tmp <- sapply(seq(0.55,0.95,0.05), function(prob) {
     return(length(stability.selection(stab1, piThreshold=prob, rhoError=10)$selected.predictors))
})
plot(seq(0.55,0.95,0.05), tmp)

# effect of restricting hyper-paramter grid
tmp <- sapply(seq(1,100,1), function(r) {
     return(length(stability.selection(stab1, piThreshold=0.75, rhoError=r)$selected.predictors))
})
plot(seq(1,100,1), tmp)

# number of selected variables
length(sample1$sel)
tmp
