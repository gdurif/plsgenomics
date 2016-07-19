
###### testing rirls.spls
     
# sources
source("pkg/R/sample.bin.R")
source("pkg/R/wirrls.R")
source("pkg/R/rpls.R")

# sample
n = 100
p = 50
kstar = 12
lstar = 12
beta.min = 0.5
beta.max = 1
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

##### test without Xtest

model1 = rpls(Ytrain=Y, Xtrain=X[,1], Lambda=2, ncomp=2, Xtest=NULL, NbIterMax=50)

str(model1)

sum(model1$Ytrain!=model1$hatY)
plot(model1$X.score[,1:2])
plot(model1$Coefficients)

sum(sort(model1$A) %in% sort(sample1$sel))/sample1$p0

##### test with Xtest

model2 = rpls(Ytrain=Y[1:80,], Xtrain=X[1:80,1], Lambda=2, ncomp=2, Xtest=X[81:100,1], NbIterMax=50)

str(model2)

sum(Y[81:100,]!=model2$hatYtest)
sum(model2$Ytrain!=model2$hatY)
plot(model2$X.score[,1:2])
plot(model2$Coefficients)

sort(model2$A)
sample1$sel
