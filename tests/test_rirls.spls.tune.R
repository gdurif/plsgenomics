###### testing rirls.spls.tune

# sources
source("pkg/R/sample.bin.R")
source("pkg/R/spls.adapt.R")
source("pkg/R/spls.adapt.aux.R")
source("pkg/R/ust.adapt.R")
source("pkg/R/ust.R")
source("pkg/R/wpls.R")

home = Sys.getenv("HOME")
source("pkg/R/wirrls.R")

library(parallel)
library(boot)

source("pkg/R/rirls.spls.R")
source("pkg/R/rirls.spls.tune.R")

# sample
n = 100
p = 50
kstar = 10
lstar = 3
beta.min = 0.5
beta.max = 1
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

##### tuning parameter

cv1 = rirls.spls.tune(X=X, Y=Y, lambda.ridge.range=1:2, lambda.l1.range=seq(0.05,0.95,by=0.3), ncomp.range=1:2, adapt=FALSE, maxIter=100, svd.decompose=TRUE, return.grid=TRUE, ncores=4, nfolds=10, center.X=FALSE, scale.X=FALSE, weighted.center=FALSE)

str(cv1)
