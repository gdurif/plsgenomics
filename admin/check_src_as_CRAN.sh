#!/bin/bash

R CMD check --as-cran plsgenomics_$(cat ../Version).tar.gz
