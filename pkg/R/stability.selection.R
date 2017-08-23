### stability.selection.R  (2016-12)
###
###
### Copyright 2016-12 Ghislain DURIF
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

#' @title
#' Stability selection procedure to select covariates for the LOGIT-SPLS 
#' and multinomial-SPLS methods
#' @aliases stability.selection
#' 
#' @description 
#' The function \code{stability.selection} returns the list of selected 
#' covariates, when following the stability selection procedure described in 
#' Durif et al. (2017). In particular, it selects covariates that are selected 
#' by most of the logit-SPLS or multinomial-SPLS models when exploring the  
#' grid of hyper-parameter candidate values.
#' 
#' 
#' @details
#' The procedure is described in Durif et al. (2017). The stability selection 
#' procedure can be summarize as follow (c.f. Meinshausen and Buhlmann, 2010).
#' 
#' (i) For each candidate values \code{(ncomp, lambda.l1, lambda.ridge)} of 
#' hyper-parameters, a logit-SPLS is trained on \code{nresamp} resamplings 
#' of the data. Then, for each triplet \code{(ncomp, lambda.l1, lambda.ridge)}, 
#' the probability that a covariate (i.e. a column in \code{X}) is selected is 
#' computed among the resamplings.
#' 
#' The estimated probabilities can be visualized as a heatmap with the 
#' function \code{\link{stability.selection.heatmap}}.
#' 
#' (ii) Eventually, the set of "stable selected" variables corresponds to the 
#' set of covariates that were selected by most of the training among the 
#' grid of hyper-parameters candidate values, based on a threshold probability
#' \code{piThreshold} and a restriction of the grid of hyper-parameters based 
#' on \code{rhoError} (c.f. Durif et al., 2017, for details).
#' 
#' This function achieves the second step (ii) of the stability selection 
#' procedure. The first step (i) is achieved by the functions
#' \code{\link{logit.spls.stab}} or \code{\link{multinom.spls.stab}}.
#' 
#' @param stab.out the output of the functions \code{\link{logit.spls.stab}}
#' or \code{\link{multinom.spls.stab}}.
#' @param piThreshold a value in (0,1], corresponding to the threshold 
#' probability used to select covariate (c.f. Durif et al., 2017).
#' @param rhoError a positive value used to restrict the grid of 
#' hyper-parameter candidate values (c.f. Durif et al., 2017).
#' 
#' @return the list of the name of covariates that are selected.
#' 
#' @author
#' Ghislain Durif (\url{http://thoth.inrialpes.fr/people/gdurif/}).
#' 
#' @seealso \code{\link{logit.spls.stab}}, \code{\link{multinom.spls.stab}}, 
#' \code{\link{stability.selection.heatmap}}
#' 
#' @examples
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 100
#' sample1 <- sample.bin(n=n, p=p, kstar=10, lstar=2, 
#'                       beta.min=0.25, beta.max=0.75, mean.H=0.2, 
#'                       sigma.H=10, sigma.F=5)
#' 
#' X <- sample1$X
#' Y <- sample1$Y
#' 
#' ### pertinent covariates id
#' sample1$sel
#' 
#' ### hyper-parameters values to test
#' lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
#' ncomp.range <- 1:10
#' # log-linear range between 0.01 a,d 1000 for lambda.ridge.range
#' logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
#' lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)
#' 
#' ### tuning the hyper-parameters
#' stab1 <- logit.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
#'                          lambda.l1.range=lambda.l1.range, 
#'                          ncomp.range=ncomp.range, 
#'                          adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
#'                          ncores=1, nresamp=100)
#'                        
#' str(stab1)
#' 
#' ### selected covariates
#' stability.selection(stab1, piThreshold=0.6, rhoError=10)
#' 
#' @export
stability.selection <- function(stab.out, piThreshold=0.9, rhoError=10) {
     
     q_Lambda <- stab.out$q.Lambda
     probs_lambda <- stab.out$probs.lambda
     
     p <- ncol(probs_lambda[,-c(1:3)])
     
     ## select variables
     q_LambdaMax <- sqrt((2*piThreshold-1) * p * rhoError)
     
     lambda_ok <- (1:nrow(q_Lambda))[which(q_Lambda$qLambda <= q_LambdaMax)]
     
     tmp_probs <- as.matrix(probs_lambda[lambda_ok,-c(1:3)])
     
     which_var <- apply(tmp_probs, 2, function(x) max(x) >= piThreshold)
     
     selected.variables <- colnames(probs_lambda[,-1])[which_var]
     
     return(selected.variables)
}

#' @title
#' Heatmap visualization of estimated probabilities of selection for each
#' covariate
#' @aliases stability.selection.heatmap
#' 
#' @description 
#' The function \code{stability.selection.heatmap} allows to visualize 
#' estimated probabilities to be selected for each covariate depending on the
#' value of hyper-parameters in the logit-spls or multinomial-spls model. 
#' These estimated probabilities are used in the stability selection procedure 
#' described in Durif et al. (2017).
#' 
#' 
#' @details
#' The procedure is described in Durif et al. (2017). The stability selection 
#' procedure can be summarize as follow (c.f. Meinshausen and Buhlmann, 2010).
#' 
#' (i) For each candidate values \code{(ncomp, lambda.l1, lambda.ridge)} of 
#' hyper-parameters, a logit-SPLS is trained on \code{nresamp} resamplings 
#' of the data. Then, for each triplet \code{(ncomp, lambda.l1, lambda.ridge)}, 
#' the probability that a covariate (i.e. a column in \code{X}) is selected is 
#' computed among the resamplings.
#' 
#' The estimated probabilities can be visualized as a heatmap with the 
#' function \code{\link{stability.selection.heatmap}}.
#' 
#' (ii) Eventually, the set of "stable selected" variables corresponds to the 
#' set of covariates that were selected by most of the training among the 
#' grid of hyper-parameters candidate values, based on a threshold probability
#' \code{piThreshold} and a restriction of the grid of hyper-parameters based 
#' on \code{rhoError} (c.f. Durif et al., 2017, for details).
#' 
#' This function allows to visualize probabalities estimated at the first 
#' step (i) of the stability selection by the functions
#' \code{\link{logit.spls.stab}} or \code{\link{multinom.spls.stab}}.
#' 
#' This function use the function \code{\link{matrix.heatmap}}.
#' 
#' @param stab.out the output of the functions \code{\link{logit.spls.stab}}
#' or \code{\link{multinom.spls.stab}}.
#' @param ... any argument that could be pass to the functions 
#' \code{\link[fields]{image.plot}} or \code{\link[graphics]{image}}.
#' 
#' @return No return, just plot the heatmap in the current graphic window.
#' 
#' @author
#' Ghislain Durif (\url{http://thoth.inrialpes.fr/people/gdurif/}).
#' 
#' @seealso \code{\link{logit.spls}}, \code{\link{stability.selection}}, 
#' \code{\link{stability.selection.heatmap}}
#' 
#' @examples
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 100
#' sample1 <- sample.bin(n=n, p=p, kstar=10, lstar=2, 
#'                       beta.min=0.25, beta.max=0.75, mean.H=0.2, 
#'                       sigma.H=10, sigma.F=5)
#' 
#' X <- sample1$X
#' Y <- sample1$Y
#' 
#' ### pertinent covariates id
#' sample1$sel
#' 
#' ### hyper-parameters values to test
#' lambda.l1.range <- seq(0.05,0.95,by=0.1) # between 0 and 1
#' ncomp.range <- 1:10
#' # log-linear range between 0.01 a,d 1000 for lambda.ridge.range
#' logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n))
#' lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=21), digits=3)
#' 
#' ### tuning the hyper-parameters
#' stab1 <- logit.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
#'                          lambda.l1.range=lambda.l1.range, 
#'                          ncomp.range=ncomp.range, 
#'                          adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
#'                          ncores=1, nresamp=100)
#'                        
#' str(stab1)
#' 
#' ### heatmap of estimated probabilities
#' stability.selection.heatmap(stab1)
#' 
#' @export
stability.selection.heatmap <- function(stab.out, ...) {
     matrix.heatmap(as.matrix(stab.out$probs.lambda[,-c(1:3)]),
                    xlab="covariates", ylab="selection sharpness",
                    ...)
}
