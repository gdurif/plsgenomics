###### testing spls.adapt

# sources

RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("sourceDir.R")

sourceDir("pkg/R")


## data
n = 500
p = 5000
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


## function parameters
Xtrain=X
Ytrain=Y
ncomp=2
weight.mat=NULL
Xtest=NULL
type="pls1"
center.X=TRUE
scale.X=FALSE
center.Y=TRUE
scale.Y=FALSE
weighted.center=FALSE


wpls(Xtrain, Ytrain, ncomp, weight.mat, Xtest, 
     type, 
     center.X, scale.X, center.Y, scale.Y, 
     weighted.center)
