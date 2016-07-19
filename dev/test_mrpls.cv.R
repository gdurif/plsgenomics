
###### testing rirls.spls

# sources
source("pkg/R/sample.multinom.R")
source("pkg/R/mwirrls.R")
source("pkg/R/mrpls.cv.R")
source("pkg/R/mrplsaux.R")

library(reshape2)

# sample
n = 100
p = 50
nb.class=4
kstar = 12
lstar = 12
beta.min = 0.05
beta.max = 0.1
mean.H=0
sigma.H=3
mean.F=0
sigma.F=1

sample1 = sample.multinom(n, p, nb.class, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

##### test without Xtest

cv1 = mrpls.cv(Ytrain=Y, Xtrain=X, LambdaRange=c(0.1, 1, 5, 10), ncompMax=10, NbIterMax=50, ncores=8)

str(cv1)
