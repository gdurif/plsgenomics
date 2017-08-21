###### testing rirls.spls.tune

# sources
library(parallel)
library(boot)
library(plyr)

setwd("/home/durif/source_code/plsgenomics")

source("pkg/R/sample.multinom.R")
source("pkg/R/spls.adapt.R")
source("pkg/R/spls.in.R")
source("pkg/R/spls.adapt.aux.R")
source("pkg/R/ust.adapt.R")
source("pkg/R/ust.R")
source("pkg/R/wpls.R")

source("pkg/R/wirrls.R")

source("pkg/R/rirls.spls.R")
source("pkg/R/rirls.spls.aux.R")
source("pkg/R/rirls.spls.stab.R")
source("pkg/R/stab_sel.R")
source("pkg/R/rirls.spls.tune2.R")

source("pkg/R/spls.in.R")
source("pkg/R/mwirrls.R")
source("pkg/R/m.rirls.spls.R")
source("pkg/R/m.rirls.spls.aux.R")
source("pkg/R/m.rirls.spls.tune2.R")
source("pkg/R/m.rirls.spls.stab.R")

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
Y = sample1$Y

print(table(Y))

##### tuning parameter
lambda.ridge.range=c(0.1, 1, 5, 10)
lambda.l1.range=seq(0.05,0.95,by=0.3)
ncomp=1
piThreshold=0.6
rhoError=10
adapt=FALSE
maxIter=100
svd.decompose=TRUE
return.grid=FALSE
ncores=2
nresamp=2
nfolds=5
center.X=TRUE
scale.X=FALSE
weighted.center=TRUE
seed=NULL
verbose=TRUE

time1 <- system.time( stab1 <- m.rirls.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, lambda.l1.range=lambda.l1.range, ncomp=1, 
                                                 adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                                                 ncores=5, nresamp=50, nfolds=5, 
                                                 center.X=TRUE, scale.X=FALSE, weighted.center=TRUE, 
                                                 seed=NULL, verbose=TRUE) )

str(stab1)

stab_sel(stab_out=stab1, piThreshold=0.6, rhoError=150)

sample1$sel

time1
