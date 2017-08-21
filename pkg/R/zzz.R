.onLoad <- function(libname, pkgname) {
     print("'plsganomics' package")
     print("C++ based sparse PLS routines will soon be available on the CRAN in the new 'fastPLS' package.")
     
     ## setting BLAS and OPENMP number of threads (modified in each function)
     blas_set_num_threads(1)
     omp_set_num_threads(1)
}
