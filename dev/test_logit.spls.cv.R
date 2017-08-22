###### testing logit.spls.cv

# sources
RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("env.R")

# sample
n = 100
p = 100
kstar = 10
lstar = 2
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

### hyper-parameters values to test
lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
ncomp.range <- 1:5
# log-linear range between 0.01 a,d 1000 for lambda.ridge.range
logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)

##### tuning parameter

time1 <- system.time( cv1 <- logit.spls.cv(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
                                           lambda.l1.range=lambda.l1.range, 
                                           ncomp.range=ncomp.range, adapt=FALSE, maxIter=100, svd.decompose=FALSE, return.grid=TRUE, 
                                           ncores=8, nfolds=5, nrun=1, center.X=TRUE, scale.X=FALSE, weighted.center=FALSE, seed=1) )

str(cv1)

time1
