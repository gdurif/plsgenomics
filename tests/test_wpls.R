###### testing wpls(Xtrain, Ytrain, ncomp, weight.mat=NULL, Xtest=NULL, type="pls1", center.X=TRUE, scale.X=FALSE, center.Y=TRUE, scale.Y=FALSE, weighted.center=FALSE)

# sources
source("R/sample.cont.R")
source("R/spls.adapt.R")
source("R/ust.adapt.R")
source("R/wpls.R")

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

### test without Xtest

model1 = wpls(Xtrain=X, Ytrain=Y, ncomp=12, weight.mat=NULL, Xtest=NULL, type="pls1", center.X=TRUE, scale.X=TRUE, center.Y=TRUE, scale.Y=TRUE, weighted.center=FALSE)

str(model1)

plot(model1$Ytrain, model1$hatY.nc)
plot(model1$sYtrain, model1$residuals)
plot(model1$X.score[,1:2])
plot(model1$betahat.nc)

### test with Xtest

model2 = wpls(Xtrain=X[1:80,], Ytrain=Y[1:80,], ncomp=2, weight.mat=NULL, Xtest=X[81:100,], type="pls1", center.X=TRUE, scale.X=TRUE, center.Y=TRUE, scale.Y=TRUE, weighted.center=FALSE)

str(model2)

plot(Y[81:100,], model2$hatYtest.nc)
points(-2000:2000, -2000:2000, type="l")
plot(model2$sYtrain, model2$residuals)
plot(model1$X.score[,1:2])
plot(model1$betahat.nc)

### with weighted scalar product

model3 = wpls(Xtrain=X[1:80,], Ytrain=Y[1:80,], ncomp=2, weight.mat=diag(1, nrow=80, ncol=80), Xtest=X[81:100,], type="pls1", center.X=TRUE, scale.X=TRUE, center.Y=TRUE, scale.Y=TRUE, weighted.center=TRUE)

str(model3)

plot(Y[81:100,], model3$hatYtest.nc)
points(-2000:2000, -2000:2000, type="l")
plot(model3$sYtrain, model3$residuals)
plot(model1$X.score[,1:2])
plot(model1$betahat.nc)

