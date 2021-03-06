###### testing rirls.spls

rm(list=ls())

# sources
RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("env.R")


# sample
n = 100
p = 50
nb.class = 4
kstar = 10
lstar = 2
beta.min = 0.05
beta.max = 0.1
mean.H=0
sigma.H=5
mean.F=0
sigma.F=1

sample1 = sample.multinom(n, p, nb.class, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

print(table(Y))

### test

time1 <- system.time( cv1 <- multinom.spls.cv(X=X, Y=Y, lambda.ridge.range=c(1, 5, 10), lambda.l1.range=seq(0.05,0.95,by=0.3), 
                                              ncomp.range=1:3, adapt=FALSE, maxIter=100, svd.decompose=TRUE, return.grid=TRUE, 
                                              ncores=8, nfolds=5, nrun=1, center.X=TRUE, scale.X=FALSE, weighted.center=FALSE, seed=1) )


str(cv1)

time1

cv1$cv.grid

### model

model1 = multinom.spls(Xtrain=X, Ytrain=Y, lambda.ridge=cv1$lambda.ridge.opt, 
                       lambda.l1=cv1$lambda.l1.opt, ncomp=cv1$ncomp.opt, 
                       Xtest=NULL, adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                       center.X=TRUE, scale.X=TRUE, weighted.center=TRUE)

sum(model1$Ytrain!=model1$hatY)/length(Y)

values=c(sample1$B[,1], model1$Coefficients[-1,1], model1$hatBeta[-1,1])
plot(sample1$B[,1], ylim=c(min(values), max(values)))
points(model1$Coefficients[-1,1], col="red", pch=2)
points(model1$hatBeta[-1,1], col="blue", pch=3)

values=c(sample1$B[,2], model1$Coefficients[-1,2], model1$hatBeta[-1,2])
plot(sample1$B[,2], ylim=c(min(values), max(values)))
points(model1$Coefficients[-1,2], col="red", pch=2)
points(model1$hatBeta[-1,2], col="blue", pch=3)


