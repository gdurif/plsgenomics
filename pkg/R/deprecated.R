### deprecated.R  (2014-10)
###
###    Deprecated functions (that are now renamed) from the 'plsgenomics'
###    package
###
### Copyright 2017-08 Ghislain DURIF
###
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

#' Deprecated function(s) in the 'plsgenomics' package
#' 
#' These functions are provided for compatibility with older version of
#' the 'plsgenomics' package.  They may eventually be completely
#' removed.
#' @rdname plsgenomics-deprecated
#' @name plsgenomics-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @docType package
#' @export  spls.adapt.tune
#' @aliases spls.adapt.tune
#' @section Details:
#' \tabular{rl}{
#'     \code{spls.adapt.tune} \tab is replaced by \code{\link{spls.cv}}\cr
#' }
#'  
spls.adapt.tune <- function(...) {
     .Deprecated("spls.cv", package="plsgenomics")
     spls.cv(...)
}
NULL