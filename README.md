# R-package plsgenomics

PLS Analyses for Genomics


## Author

Anne-Laure Boulesteix <boulesteix@ibe.med.uni-muenchen.de>, Ghislain Durif <gd.dev@libertymail.net>, Sophie Lambert-Lacroix <sophie.lambert-lacroix@univ-grenoble-alpes.fr>, Julie Peyre <Julie.Peyre@univ-grenoble-alpes.fr>, and Korbinian Strimmer <k.strimmer@imperial.ac.uk>.

Maintainer: Ghislain Durif <gd.dev@libertymail.net>


## Description

Routines for PLS-based genomic analyses, implementing PLS methods for classification with microarray data and prediction of transcription factor activities from combined ChIP-chip analysis. The >=1.2-1 versions include two new classification methods for microarray data: GSIM and Ridge PLS. The >=1.3 versions includes a new classification method combining variable selection and compression in logistic regression context: logit-SPLS; and an adaptive version of the sparse PLS.


## Installation

You can install the `plsgenomics` R package with the following R commands:

```R
install.packages("plsgenomics")
```

or

```R
devtools::install_git("https://gitlab.inria.fr/gdurif/pCMF", subdir="pkg")
```

To install the `devtools` package, you can run:

```R
install.packages("devtools")
```

You can also use the git repository available at <https://gitlab.inria.fr/gdurif/plsgenomics>, 
then build and install the package with Rstudio (the [project file](./plsgenomics.Rproj) 
is set accordingly) or with the R command line tools.

Once you cloned the git repository, you can also run:

```R
devtools::install("/path/to/plsgenomics/pkg") # you should edit the path
```


## Licence

The `plsgenomics` package is distributed under the GPL (>=2) licence.


## Example

Examples regarding the sparse PLS method and the sparse PLS approach for logistic regression developped in Durif et al. (2018) can be respectively found in these two scripts: [spls_example.R](./spls_example.R) and [logit_spls_example.R](./logit_spls_example.R).


## Reference

Durif, G., Modolo, L., Michaelsson, J., Mold, J.E., Lambert-Lacroix, S., Picard, F., 2018. High dimensional classification with combined adaptive sparse PLS and logistic regression. Bioinformatics 34, 485–493. https://doi.org/10.1093/bioinformatics/btx571

