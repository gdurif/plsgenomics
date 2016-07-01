###### testing mrpls.cv

rm(list=ls())

# sources
source("pkg/R/mrpls.cv.R")
source("pkg/R/mrplsaux.R")
source("pkg/R/mwirrls.R")
source("pkg/R/sample.multinom.R")

library(reshape2)

# sample
n = 100
p = 50
nb.class=3
kstar = 12
lstar = 3
beta.min = 0.05
beta.max = 0.1
mean.H=0
sigma.H=3
mean.F=0
sigma.F=1

sample1 = sample.multinom(n, p, nb.class, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = as.vector(sample1$Y)


mrpls.cv(Ytrain=Y, Xtrain=X, LambdaRange=c(0.1,1),ncompMax=3, ncores=5)
