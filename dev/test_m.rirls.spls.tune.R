###### testing rirls.spls

rm(list=ls())

# sources
# source("pkg/R/spls.in.R")
# source("pkg/R/mwirrls.R")
# source("pkg/R/m.rirls.spls.R")
# source("pkg/R/m.rirls.spls.aux.R")
# source("pkg/R/m.rirls.spls.tune.R")
# source("pkg/R/ust.adapt.R")
# source("pkg/R/ust.R")
# source("pkg/R/wpls.R")
# source("pkg/R/sample.multinom.R")
# 
# # library
# library(parallel)
# library(MASS)
library(plsgenomics)


# sample
n = 100
p = 50
nb.class=3
kstar = 12
lstar = 3
beta.min = 0.5
beta.max = 1
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.multinom(n, p, nb.class, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

### test

cv1 = m.rirls.spls.tune(X=X, Y=Y, lambda.ridge.range=c(1,10), lambda.l1.range=seq(0.05,0.95,by=0.3), ncomp.range=2:3, adapt=FALSE, maxIter=100, svd.decompose=TRUE, return.grid=TRUE, ncores=4, nfolds=5, center.X=TRUE, scale.X=TRUE, weighted.center=TRUE, nrun=5, verbose=TRUE)

str(cv1)
