###### testing spls.cv

# sources
RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("env.R")

# sample
n = 100
p = 1000
kstar = 500
lstar = 250
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

cv1 = spls.cv(X=X, Y=Y, lambda.l1.range=seq(0.05, 0.95, by=0.1), ncomp.range=1:10, 
              weight.mat=NULL, adapt=FALSE, center.X=TRUE, center.Y=TRUE, 
              scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, 
              return.grid=TRUE, ncores=8, nfolds=10, nrun=1)

str(cv1)
