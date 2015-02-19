###### testing spls.adapt.tune

# sources
source("R/sample.cont.R")
source("R/spls.adapt.aux.R")
source("R/spls.adapt.tune.R")
source("R/ust.adapt.R")
source("R/ust.R")
source("R/wpls.R")

library(parallel)

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

##### model tuning

cv1 = spls.adapt.tune(X=X, Y=Y, lambda.l1.range=seq(0.05, 0.95, by=0.3), ncomp.range=1:2, weight.mat=NULL, adapt=TRUE, center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, return.grid=TRUE, ncores=4, nfolds=10)

str(cv1)
