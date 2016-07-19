###### testing rirls.spls.tune

# sources
library(parallel)
library(boot)

setwd("/home/durif/source_code/plsgenomics")

source("pkg/R/sample.bin.R")
source("pkg/R/spls.adapt.R")
source("pkg/R/spls.in.R")
source("pkg/R/spls.adapt.aux.R")
source("pkg/R/ust.adapt.R")
source("pkg/R/ust.R")
source("pkg/R/wpls.R")

source("pkg/R/wirrls.R")

source("pkg/R/rirls.spls.R")
source("pkg/R/rirls.spls.aux.R")
source("pkg/R/rirls.spls.tune.R")
source("pkg/R/rirls.spls.tune2.R")

# sample
n = 50
p = 50
kstar = 10
lstar = 10
beta.min = 5
beta.max = 10
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

print(table(Y))

##### tuning parameter

time1 <- system.time( cv1 <- rirls.spls.tune(X=X, Y=Y, lambda.ridge.range=c(0.1, 1, 5, 10), lambda.l1.range=seq(0.05,0.95,by=0.1), 
                                             ncomp.range=1:4, adapt=FALSE, maxIter=100, svd.decompose=FALSE, return.grid=TRUE, 
                                             ncores=8, nfolds=5, nrun=1, center.X=TRUE, scale.X=FALSE, weighted.center=FALSE, seed=1) )


time2 <- system.time( cv2 <- rirls.spls.tune2(X=X, Y=Y, lambda.ridge.range=c(0.1, 1, 5, 10), lambda.l1.range=seq(0.05,0.95,by=0.1), 
                                              ncomp.range=1:4, adapt=FALSE, maxIter=100, svd.decompose=FALSE, return.grid=TRUE, 
                                              ncores=8, nfolds=5, nrun=1, center.X=TRUE, scale.X=FALSE, weighted.center=FALSE, seed=1) )

str(cv1)
str(cv2)

time1
time2
