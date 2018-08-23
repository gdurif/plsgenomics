# how to use logit-spls

library(plsgenomics)

### generate some data sample
## X: matrix n x p with observations in rows and covariates in column
## Y: vector of length n with binary values (0-1)
n = 200
p = 100
kstar = 10
lstar = 2
beta.min = 5
beta.max = 10
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

print(table(Y))

### splitting between learning and testing set (70/30%)
index.train <- sort(sample(1:n, size=round(0.7*n)))
index.test <- (1:n)[-index.train]

Xtrain <- X[index.train,]
Ytrain <- Y[index.train,]

Xtest <- X[index.test,]
Ytest <- Y[index.test,]


######## tunning hyper-parameters with cross-validation
### hyper-parameters values to test
lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
ncomp.range <- 1:5
# log-linear range between 0.01 a,d 1000 for lambda.ridge.range
logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)

### tuning the hyper-parameters 
## (/!\ using 8 cores)
cv1 <- logit.spls.cv(X=Xtrain, Y=Ytrain, lambda.ridge.range=lambda.ridge.range, 
                     lambda.l1.range=lambda.l1.range, 
                     ncomp.range=ncomp.range, adapt=TRUE, maxIter=100, svd.decompose=TRUE, return.grid=TRUE, 
                     ncores=8, nfolds=5, nrun=1, center.X=TRUE, scale.X=FALSE, weighted.center=TRUE)

str(cv1)


### otpimal values
lambda.l1 <- cv1$lambda.l1.opt
lambda.ridge <- cv1$lambda.ridge.opt
ncomp <- cv1$ncomp.opt

######## fitting the model, and predicting new observations
model1 <- logit.spls(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, 
                     lambda.ridge=lambda.ridge, lambda.l1=lambda.l1, ncomp=ncomp, 
                     adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                     center.X=TRUE, scale.X=FALSE, weighted.center=TRUE)
str(model1)

### estimation error
sum(Ytrain!=model1$hatY) / length(Ytrain)
### prediction error
sum(Ytest!=model1$hatYtest) / length(Ytest)

### plotting the predictor coefficients
plot(model1$betahat.nc, xlab="variable index", ylab="coeff")

### low dimensional space representation
plot(model1$X.score[,1:2], xlab="comp1", ylab="comp2", col=Ytrain+1)

### selected variables
model1$A



######### stability selection procedure
## (/!\ using 8 cores)
stab1 <- logit.spls.stab(X=Xtrain, Y=Ytrain, lambda.ridge.range=lambda.ridge.range, 
                         lambda.l1.range=lambda.l1.range, 
                         ncomp.range=ncomp.range, 
                         adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                         ncores=8, nresamp=50, 
                         center.X=TRUE, scale.X=FALSE, weighted.center=TRUE, 
                         seed=NULL, verbose=TRUE)

str(stab1)

### heatmap of estimated probabilities
stability.selection.heatmap(stab1)

### selected covariates
tmp <- stability.selection(stab1, piThreshold=0.75, rhoError=10)
tmp

### "true" pertinent covariates
sample1$sel
