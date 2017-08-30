###### testing logit.spls.stab

# sources
RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("env.R")

# ## debugging
# # sample
# n = 50
# p = 50
# kstar = 10
# lstar = 2
# beta.min = 5
# beta.max = 10
# mean.H=0
# sigma.H=10
# mean.F=0
# sigma.F=5
# 
# sample1 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)
# 
# X = sample1$X
# Y = sample1$Y
# 
# print(table(Y))
# 
# ##### tuning parameter
# lambda.ridge.range=c(0.1, 1, 5, 10)
# lambda.l1.range=seq(0.05,0.95,by=0.1)
# ncomp.range=c(1,2)
# piThreshold=0.6
# rhoError=10
# adapt=FALSE
# maxIter=100
# svd.decompose=TRUE
# return.grid=FALSE
# ncores=2
# nresamp=40
# nfolds=5
# center.X=TRUE
# scale.X=FALSE
# weighted.center=TRUE
# seed=NULL
# verbose=TRUE
# 
# time1 <- system.time( stab1 <- logit.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range,
#                                                lambda.l1.range=lambda.l1.range, ncomp.range=1:5,
#                                                adapt=TRUE, maxIter=100, svd.decompose=TRUE,
#                                                ncores=8, nresamp=50,
#                                                center.X=TRUE, scale.X=FALSE, weighted.center=TRUE,
#                                                seed=NULL, verbose=TRUE) )
# 
# str(stab1)
# stab.sel.heatmap(stab1)
# 
# stab_sel(stab_out=stab1, piThreshold=0.6, rhoError=15)
# 
# sample1$sel
# 
# time1


### full example
### generating data
n <- 100
p <- 50
sample1 <- sample.bin(n=n, p=p, kstar=10, lstar=2,
                      beta.min=0.25, beta.max=0.75, mean.H=0.2,
                      sigma.H=10, sigma.F=5)

X <- sample1$X
Y <- sample1$Y

sample1$sel

### hyper-parameters values to test
lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
ncomp.range <- 1:5
# log-linear range between 0.01 a,d 1000 for lambda.ridge.range
logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)

### tuning the hyper-parameters
time1 <- system.time( stab1 <- logit.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, lambda.l1.range=lambda.l1.range, 
                                               ncomp.range=ncomp.range, 
                                               adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                                               ncores=12, nresamp=50, 
                                               center.X=TRUE, scale.X=FALSE, weighted.center=TRUE, 
                                               seed=NULL, verbose=TRUE))

str(stab1)

time1

### heatmap of estimated probabilities
stability.selection.heatmap(stab1)

### selected covariates
tmp <- stability.selection(stab1, piThreshold=0.75, rhoError=10)
tmp
