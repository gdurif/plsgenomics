
###### testing rirls.spls

# sources
source("pkg/R/sample.bin.R")
source("pkg/R/wirrls.R")
source("pkg/R/rpls.cv.R")
source("pkg/R/rplsaux.R")

library(reshape2)

# sample
n = 300
p = 500
kstar = 12
lstar = 1
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

cv1 = rpls.cv(Ytrain=Y, Xtrain=X, LambdaRange=c(1,5,10,100), ncompMax=4, NbIterMax=50, ncores=4)

str(cv1)
