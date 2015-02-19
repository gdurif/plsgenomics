##### testing package installation


system("R CMD INSTALL --library=/home/durif/source_code/tests plsgenomics_1.3.tar.gz")

library(plsgenomics, lib.loc="/home/durif/source_code/tests")
