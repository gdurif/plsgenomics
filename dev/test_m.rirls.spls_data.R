###### test mspls on specific data to fix a bug
rm(list=ls())

# sources
source("pkg/R/spls.in.R")
source("pkg/R/mwirrls.R")
source("pkg/R/m.rirls.spls.R")
source("pkg/R/m.rirls.spls.aux.R")
source("pkg/R/m.rirls.spls.tune.R")
source("pkg/R/ust.adapt.R")
source("pkg/R/ust.R")
source("pkg/R/wpls.R")
source("pkg/R/sample.multinom.R")

# library
library(parallel)
library(MASS)

## files
wd = "/home/durif/script_code/2015/singleCell/"

dataDir = "/storage/durif/data/2015/singleTcell/v1/"


load(paste0(dataDir, "mspls_test.RData"))
load(paste0(dataDir, "tmp_mspls.RData"))
r_select_a <- data_infos$day == "D15" & is.na(data_infos$phenotype)
r_select_b <- data_infos$day == "D15" & !is.na(data_infos$phenotype) 
print(dim(as.matrix(data_genes_norm_batch[r_select_b,c_select])))
print(dim(as.matrix(data_genes_norm_batch[r_select_a,c_select])))
print(length(tmp_mspls$parameters$is_group))
tmp_mspls$parameters$is_group
model <- m.rirls.spls(Xtrain=as.matrix(data_genes_norm_batch[r_select_b,c_select]),
                      Ytrain=tmp_mspls$parameters$is_group,
                      lambda.ridge=tmp_mspls$parameters$lambda.ridge,
                      lambda.l1=tmp_mspls$parameters$lambda.l1,
                      ncomp=tmp_mspls$parameters$ncomp,
                      Xtest=as.matrix(data_genes_norm_batch[r_select_a,c_select]),
                      adapt=TRUE,
                      maxIter=100,
                      svd.decompose=FALSE,
                      center.X=TRUE,
                      scale.X=TRUE,
                      weighted.center=TRUE)


str(as.matrix(data_genes_norm_batch[r_select_b,c_select]))
