###### testing rirls.spls

rm(list=ls())

# sources
source("pkg/R/sample.multinom.R")
source("pkg/R/mwirrls.R")
source("pkg/R/mrpls.R")

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

model1 = mrpls(Ytrain=Y, Xtrain=X, Lambda=2, ncomp=2, Xtest=NULL, NbIterMax=50)

str(model1)

sum(model1$Ytrain!=model1$hatY)
plot(model1$X.score[,1:2])
plot(model1$Coefficients)

##### test with Xtest

model2 = mrpls(Ytrain=Y[1:80,], Xtrain=X[1:80,1], Lambda=2, ncomp=2, Xtest=X[81:100,1], NbIterMax=50)

str(model2)

sum(Y[81:100,]!=model2$hatYtest)
sum(model2$Ytrain!=model2$hatY)
plot(model2$X.score[,1:2])
plot(model2$Coefficients)

sort(model2$A)
sample1$sel



### tests on components

G=3
ntrain = n

X.score <- lapply(1:G, function(g) {
     res <- model1$Xtrain[ (0:(ntrain-1)) * G + g, (g-1)*p + (1:p)] %*% model1$X.weight[ (g-1)*p + (1:p), ]
     return(res)
})

plot(X.score[[1]], col=(Y==1)+1)
plot(X.score[[2]], col=(Y==2)+1)
plot(X.score[[3]], col=(Y==3)+1)

X.weight <- lapply(1:G, function(g) {
     return(resSPLS$X.weight[ (g-1)*p + (1:p), ])
})


X.score2 <- lapply(1:G, function(g) {
     return(model1$X.score[ (0:(ntrain-1)) * G + g, ])
})
plot(X.score2[[1]], col=(Y==1)+1)
plot(X.score2[[2]], col=(Y==2)+1)
plot(X.score2[[3]], col=(Y==3)+1)
