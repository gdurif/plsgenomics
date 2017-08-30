###### testing rirls.spls

rm(list=ls())

# sources
RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("env.R")


# sample
n = 100
p = 50
nb.class=3
kstar = 20
lstar = 2
beta.min = 0.25
beta.max = 0.75
mean.H=0.2
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.multinom(n, p, nb.class, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

##### test without Xtest

model1 = multinom.spls(Xtrain=X, Ytrain=Y, lambda.ridge=1, lambda.l1=0.5, ncomp=2, Xtest=NULL, adapt=TRUE, maxIter=100, svd.decompose=TRUE, center.X=TRUE, scale.X=TRUE, weighted.center=TRUE)

str(model1)

sum(model1$Ytrain!=model1$hatY)/length(Y)

cbind(model1$proba, Y+1)

plot(model1$Coefficients[,1])
points(model1$Coefficients[,2], col="red")

plot(model1$hatBeta[,1])
points(model1$hatBeta[,2], col="red")

plot(sample1$B[,1])
points(sample1$B[,2], col="red")

plot(model1$Coefficients[,1], model1$hatBeta[,1])
abline(a=0,b=1)
points(model1$Coefficients[,2], model1$hatBeta[,2], col="red")

plot(sample1$B[,1], model1$hatBeta[-1,1])
abline(a=0,b=1)
points(sample1$B[,1], model1$Coefficients[-1,1], col="red")

plot(sample1$B[,1])
points(model1$Coefficients[-1,1], col="red", pch=2)
points(model1$hatBeta[-1,1], col="blue", pch=3)

print(model1$A)
print(sample1$sel)

plot(model1$X.score[[1]][,1], model1$X.score[[2]][,1], col=(Y+1))
plot(model1$X.score.full[,1], model1$X.score.full[,2], col=(Y+1))


plot(model1$X.weight[[1]][,1], model1$X.weight[[2]][,1])


##### test with Xtest

Xtrain = X[1:80,]
Xtest = X[81:100,]
Ytrain = Y[1:80,]
Ytest = Y[81:100,]

model2 = m.rirls.spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.ridge=0.0001, lambda.l1=0.3, ncomp=2, Xtest=Xtest, adapt=TRUE, maxIter=100, svd.decompose=TRUE, center.X=TRUE, scale.X=TRUE, weighted.center=FALSE)
str(model2)

sum(Ytest!=model2$hatYtest)
sum(Ytrain!=model2$hatY)
