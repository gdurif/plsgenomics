###### testing spls.stab

# sources
RDIR <- system("git rev-parse --show-toplevel", intern=TRUE)
setwd(RDIR)
source("env.R")

# sample
n = 100
p = 100
kstar = 10
lstar = 2
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

sample1$sel

##### model tuning

time1 <- system.time( stab1 <- spls.stab(X=X, Y=Y, lambda.l1.range=seq(0.05, 0.95, by=0.1), ncomp.range=1:10, 
                                         weight.mat=NULL, adapt=FALSE, center.X=TRUE, center.Y=TRUE, 
                                         scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, 
                                         ncores=8, nresamp=20))

str(stab1)

time1

### heatmap of estimated probabilities
stability.selection.heatmap(stab1)

### selected covariates
stability.selection(stab1, piThreshold=0.75, rhoError=10)
