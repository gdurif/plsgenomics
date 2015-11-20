###### testing rirls.spls.tune

# sources
library(parallel)
library(boot)

setwd("/home/durif/source_code/plsgenomics")

source("pkg/R/sample.bin.R")
source("pkg/R/spls.adapt.R")
source("pkg/R/spls.adapt.aux.R")
source("pkg/R/ust.adapt.R")
source("pkg/R/ust.R")
source("pkg/R/wpls.R")

source("pkg/R/wirrls.R")

source("pkg/R/rirls.spls.R")
source("pkg/R/rirls.spls.aux.R")
source("pkg/R/rirls.spls.tune.R")

# sample
n = 40
p = 100
kstar = 10
lstar = 3
beta.min = 5
beta.max = 10
mean.H=5
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

print(table(Y))

##### tuning parameter

cv1 = rirls.spls.tune(X=X, Y=Y, lambda.ridge.range=c(1,10), lambda.l1.range=seq(0.05,0.95,by=0.3), ncomp.range=2:3, adapt=FALSE, maxIter=100, svd.decompose=TRUE, return.grid=TRUE, ncores=4, nfolds=5, center.X=FALSE, scale.X=FALSE, weighted.center=FALSE, nrun=5)

str(cv1)
