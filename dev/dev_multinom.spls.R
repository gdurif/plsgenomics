#### dev m.rirls.spls

rm(list=ls())


## bloc definition (n <p)
Xtrain = matrix(1:10, ncol=2)
Xtest = matrix(1:4, ncol=2)
ntrain = nrow(Xtrain)
ntest=nrow(Xtest)
Z = cbind(rep(1,ntrain), Xtrain)
p <- 2
r <- p #min(p, ntrain)

G=3

Zbloc <- matrix(0,nrow=ntrain*G,ncol=G*(r+1))

if (!is.null(Xtest)) {
     Zt <- cbind(rep(1,ntest),Xtest)
     Ztestbloc <- matrix(0,nrow=ntest*G,ncol=G*(r+1))
}

for (g in 1:G) {
     row <- (0:(ntrain-1))*G+g
     col <- (r+1)*(g-1)+1:(r+1)
     Zbloc[row,col] <- Z
     if (!is.null(Xtest)) {
          row <- (0:(ntest-1))*G+g
          Ztestbloc[row,col] <- Zt
     }
}

## bloc definition (n > p)
Xtrain = matrix(1:12, ncol=4)
Xtest = matrix(1:8, ncol=4)
ntrain = nrow(Xtrain)
ntest=nrow(Xtest)
Z = cbind(rep(1,ntrain), Xtrain)
p <- 4
r <- min(p, ntrain)

G=3

Zbloc <- matrix(0,nrow=ntrain*G,ncol=G*(p+1))

if (!is.null(Xtest)) {
     Zt <- cbind(rep(1,ntest),Xtest)
     Ztestbloc <- matrix(0,nrow=ntest*G,ncol=G*(p+1))
}

for (g in 1:G) {
     row <- (0:(ntrain-1))*G+g
     col <- (p+1)*(g-1)+1:(p+1)
     Zbloc[row,col] <- Z
     if (!is.null(Xtest)) {
          row <- (0:(ntest-1))*G+g
          Ztestbloc[row,col] <- Zt
     }
}       



## dev in m.rirls.spls
Xtrain = matrix(1:15, ncol=3)
ntrain=5
Ytrain=c(0:2, 0:1)
Xtest=NULL
lambda.ridge=5
lambda.l1=0.5
NbIterMax=20
ncomp=2
adapt=TRUE
maxIter=50
getwd()
source("pkg/R/spls.in.R")
source("pkg/R/mwirrls.R")
source("pkg/R/ust.adapt.R")
source("pkg/R/ust.R")
source("pkg/R/wpls.R")
source("pkg/R/sample.multinom.R")

# sample
ntrain = 100
p = 50
nb.class=3
kstar = 12
lstar = 3
beta.min = 0.5
beta.max = 1
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.multinom(n, p, nb.class, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

Xtrain = sample1$X
Ytrain = sample1$Y
