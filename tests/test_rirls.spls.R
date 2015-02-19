###### testing rirls.spls

# sources
source("pkg/R/sample.bin.R")
source("pkg/R/spls.adapt.R")
source("pkg/R/spls.adapt.aux.R")
source("pkg/R/ust.adapt.R")
source("pkg/R/ust.R")
source("pkg/R/wpls.R")

home = Sys.getenv("HOME")
source(paste0(home, "/source_code/plsgenomics/pkg/R/wirrls.R"))

library(parallel)
library(boot)

source("pkg/R/rirls.spls.R")

# sample
n = 100
p = 1000
kstar = 12
lstar = 3
beta.min = 0.5
beta.max = 1
mean.H=0
sigma.H=10
mean.F=0
sigma.F=5

sample1 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

X = sample1$X
Y = sample1$Y

##### test without Xtest

model1 = rirls.spls(Xtrain=X, Ytrain=Y, lambda.ridge=2, lambda.l1=0.5, ncomp=2, Xtest=NULL, adapt=TRUE, maxIter=100, svd.decompose=TRUE)

str(model1)

sum(model1$Ytrain!=model1$hatY)
plot(model1$X.score[,1:2])
plot(model1$Coefficients)

sum(sort(model1$A) %in% sort(sample1$sel))/sample1$p0

##### test with Xtest

model2 = rirls.spls(Xtrain=X[1:80,], Ytrain=Y[1:80,], lambda.ridge=2, lambda.l1=0.5, ncomp=2, Xtest=X[81:100,], adapt=TRUE, maxIter=100, svd.decompose=TRUE)

str(model2)

sum(Y[81:100,]!=model2$hatYtest)
sum(model2$Ytrain!=model2$hatY)
plot(model2$X.score[,1:2])
plot(model2$Coefficients)


###########################################################
##### testing aux version
###########################################################

source("pkg/R/rirls.spls.aux.R")

Xtrain <- X[1:80,]
Ytrain <- as.matrix(Y[1:80,])

ntrain <- nrow(Xtrain)

Xtest <- X[81:100,]
Ytest <- as.matrix(Y[81:100,])

ntest <- nrow(Xtest)

r <- min(p, ntrain)
DeletedCol <- NULL

svd.decompose=TRUE

### Standardize the Xtrain matrix
# standard deviation (biased one) of Xtrain
sigma2train <- apply(Xtrain, 2, var) * (ntrain-1)/(ntrain)

# mean of Xtrain
meanXtrain <- apply(Xtrain,2,mean)

# center and scale Xtrain
sXtrain <- scale(Xtrain, center=meanXtrain, scale=sqrt(sigma2train))

sXtrain.nosvd = sXtrain # keep in memory if svd decomposition

# Compute the svd when necessary -> case p > ntrain (high dim)
if ((p > ntrain) && (svd.decompose)) {
	# svd de sXtrain
	svd.sXtrain <- svd(t(sXtrain))
	# number of singular value non null
	r <- length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
	V <- svd.sXtrain$u[,1:r]
	D <- diag(c(svd.sXtrain$d[1:r]))
	U <- svd.sXtrain$v[,1:r]
	sXtrain <- U %*% D
	rm(D)
	rm(U)
	rm(svd.sXtrain)
}

# center and scale Xtest	
meanXtest <- apply(Xtest,2,mean)
sigma2test <- apply(Xtest,2,var)

sXtest <- scale(Xtest, center=meanXtrain, scale=sqrt(sigma2train))

sXtest.nosvd <- sXtest # keep in memory if svd decomposition

# if svd decomposition
if ((p > ntrain) && (svd.decompose)) {
	sXtest <- sXtest%*%V
}


model3 = rirls.spls.aux(sXtrain=sXtrain, sXtrain.nosvd=sXtrain.nosvd, Ytrain=Ytrain, lambda.ridge=2, lambda.l1=0.5, ncomp=2, sXtest=sXtest, sXtest.nosvd=sXtest.nosvd, adapt=TRUE, maxIter=100, svd.decompose=svd.decompose, meanXtrain=meanXtrain, sigma2train=sigma2train)

str(model3)

sum(model2$hatYtest!=model3$hatYtest)
