###### testing spls.adapt

# sources
library(boot)

RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("sourceDir.R")

sourceDir("pkg/R")

# sample
n = 100
p = 1000
kstar = 12
lstar = 3
beta.min = 0.5
beta.max = 1
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5
sigma.E=15

sample1 = sample.cont(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F, sigma.E)

X = sample1$X
Y = sample1$Y


print("##### test without Xtest")

model1 = spls(Xtrain=X, Ytrain=Y, lambda.l1=0.5, ncomp=2, weight.mat=NULL, 
              Xtest=NULL, adapt=TRUE, 
              center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE)

str(model1)

# plot(model1$Ytrain, model1$hatY.nc)
# plot(model1$sYtrain, model1$residuals)
# plot(model1$X.score[,1:2])
# plot(model1$betahat.nc)


print("##### test with Xtest")

model2 = spls(Xtrain=X[1:80,], Ytrain=Y[1:80,], lambda.l1=0.5, ncomp=2, 
              weight.mat=NULL, Xtest=X[81:100,], adapt=TRUE, 
              center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE)

str(model2)

# plot(Y[81:100,], model2$hatYtest.nc)
# points(-2000:2000, -2000:2000, type="l")
# plot(model2$sYtrain, model2$residuals)
# plot(model2$X.score[,1:2])
# plot(model2$betahat.nc)

print("##### with weighted scalar product")

model3 = spls(Xtrain=X[1:80,], Ytrain=Y[1:80,], lambda.l1=0.5, ncomp=2, weight.mat=diag(1, nrow=80, ncol=80), Xtest=X[81:100,], adapt=TRUE, center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE)

str(model3)

# plot(Y[81:100,], model3$hatYtest.nc)
# points(-2000:2000, -2000:2000, type="l")
# plot(model3$sYtrain, model3$residuals)
# plot(model3$X.score[,1:2])
# plot(model3$betahat.nc)


##### test of aux version

# meanXtrain=apply(X[1:80,],2,mean)
# meanYtrain=apply(as.matrix(Y[1:80,]),2,mean)
# sigmaXtrain=apply(X[1:80,],2,sd)
# sigmaYtrain=apply(as.matrix(Y[1:80,]),2,sd)
# 
# sXtrain = scale(X[1:80,], center=meanXtrain, scale=sigmaXtrain)
# sXtest = scale(X[81:100,], center=meanXtrain, scale=sigmaXtrain)
# 
# sYtrain = scale(as.matrix(Y[1:80,]), center=meanYtrain, scale=sigmaYtrain)
# 
# 
# model4 = spls.adapt.aux(Xtrain=X[1:80,], sXtrain=sXtrain, Ytrain=Y[1:80,], sYtrain=sYtrain, lambda.l1=0.5, ncomp=2, weight.mat=NULL, Xtest=X[81:100,], sXtest=sXtest, adapt=TRUE, meanXtrain=apply(X[1:80,],2,mean), meanYtrain=apply(as.matrix(Y[1:80,]),2,mean), sigmaXtrain=apply(X[1:80,],2,sd), sigmaYtrain=apply(as.matrix(Y[1:80,]),2,sd), center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE)
# 
# str(model4)
# 
# plot(model3$hatYtest, model4$hatYtest)
# points(-2000:2000, -2000:2000, type="l")
# 
# plot(model3$betahat, model4$betahat)
# points(-2000:2000, -2000:2000, type="l")
