##### testing package installation

## clear environment
rm(list=ls())

## recreating the MD5 file
source("computeMD5sum.R")

## remove previous installation test if exists
system("if [ -d /home/durif/source_code/tests/plsgenomics ]; then rm -r /home/durif/source_code/tests/plsgenomics; fi;")

## recreating the package source archive
system("tar -zcvf plsgenomics_1.3.tar.gz pkg/")

## installing this new version of the package
system("R CMD INSTALL --library=/home/durif/source_code/tests plsgenomics_1.3.tar.gz")

## loading this new version of the package
library(plsgenomics, lib.loc="/home/durif/source_code/tests")
