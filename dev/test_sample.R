##### testing functions generating sample

# sources
source("R/sample.cont.R")
source("R/sample.bin.R")
library(boot)


# continuous response

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
sigma.E=15

sample1 = sample.cont(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F, sigma.E)

str(sample1)



# binary response

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

sample2 = sample.bin(n, p, kstar, lstar, beta.min, beta.max, mean.H, sigma.H, mean.F, sigma.F)

str(sample2)
