### spls.cv.R  (2014-10)
###
###    Tuning parameters (ncomp, lambda.l1) for adaptive sparse PLS regression 
###    for continuous response, by K-fold cross-validation
###
### Copyright 2014-10 Ghislain DURIF
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
#' Tuning parameters (ncomp, lambda.l1) for Adaptive Sparse PLS regression for 
#' continuous response, by K-fold cross-validation
#' @aliases spls.cv
#' 
#' @description 
#' The function \code{spls.cv} tuns the hyper-parameter values used 
#' in the \code{spls} procedure, by minimizing the mean squared error 
#' of prediction over the hyper-parameter grid, 
#' using Durif et al. (2017) adaptive SPLS algorithm.
#' 
#' @details 
#' The columns of the data matrices \code{Xtrain} and \code{Xtest} may not 
#' be standardized, since standardizing can be performed by the function 
#' \code{spls.cv} as a preliminary step before the algorithm is run.
#' 
#' The procedure is described in Durif et al. (2015). The K-fold 
#' cross-validation can be summarize as follow: the train set is partitioned 
#' into K folds, for each value of hyper-parameters the model is fit K times, 
#' using each fold to compute the prediction error rate, and fitting the 
#' model on the remaining observations. The cross-validation procedure returns 
#' the optimal hyper-parameters values, meaning the one that minimize 
#' the mean squared error of prediction averaged over all the folds.
#' 
#' This procedures uses the \code{mclapply} from the \code{parallel} package, 
#' available on GNU/Linux and MacOS. Users of Microsoft Windows can refer to 
#' the README file in the source to be able to use a mclapply type function.
#' 
#' @param X a (n x p) data matrix of predictors. \code{X} must be a matrix. 
#' Each row corresponds to an observation and each column to a 
#' predictor variable.
#' @param Y a (n) vector of (continuous) responses. \code{Y} must be a 
#' vector or a one column matrix. It contains the response variable for 
#' each observation.
#' @param lambda.l1.range a vecor of positive real values, in [0,1]. 
#' \code{lambda.l1} is the sparse penalty parameter for the dimension 
#' reduction step by sparse PLS (see details), the optimal values will be 
#' chosen among \code{lambda.l1.range}.
#' @param ncomp.range a vector of positive integers. \code{ncomp} is the 
#' number of PLS components. The optimal values will be chosen 
#' among \code{ncomp.range}.
#' @param weight.mat a (ntrain x ntrain) matrix used to weight the l2 metric 
#' in the observation space, it can be the covariance inverse of the Ytrain 
#' observations in a heteroskedastic context. If NULL, the l2 metric is the 
#' standard one, corresponding to homoskedastic model (\code{weight.mat} is the 
#' identity matrix).
#' @param adapt a boolean value, indicating whether the sparse PLS selection 
#' step sould be adaptive or not (see details).
#' @param center.X a boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be centered or not.
#' @param scale.X}{aa boolean value indicating whether the data matrices 
#' \code{Xtrain} and \code{Xtest} (if provided) should be scaled or not 
#' (\code{scale.X=TRUE} implies \code{center.X=TRUE}).
#' @param center.Y a boolean value indicating whether the response values 
#' \code{Ytrain} set should be centered or not.
#' @param scale.Y a boolean value indicating whether the response values 
#' \code{Ytrain} should be scaled or not (\code{scale.Y=TRUE} implies 
#' \code{center.Y=TRUE}).
#' @param weighted.center a boolean value indicating whether the centering 
#' should take into account the weighted l2 metric or not 
#' (if TRUE, it requires that weighted.mat is non NULL).
#' @param return.grid a boolean values indicating whether the grid of 
#' hyper-parameters values with corresponding mean prediction error rate over 
#' the folds should be returned or not.
#' @param ncores a positve integer, indicating the number of cores that the 
#' cross-validation is allowed to use for parallel computation (see details).
#' @param nfolds a positive integer indicating the number of folds in the 
#' K-folds cross-validation procedure, \code{nfolds=n} corresponds 
#' to the leave-one-out cross-validation, default is 10.
#' @param nrun a positive integer indicating how many times the K-folds cross-
#' validation procedure should be repeated, default is 1.
#' 
#' @return An object with the following attributes
#' \item{lambda.l1.opt}{the optimal value in \code{lambda.l1.range}.}
#' \item{ncomp.opt}{the optimal value in \code{ncomp.range}.}
#' \item{cv.grid}{the grid of hyper-parameters and corresponding prediction 
#' error rate over the folds. 
#' \code{cv.grid} is NULL if \code{return.grid} is set to FALSE.}
#' 
#' @references 
#' Durif G., Modolo L., Michaelsson J., Mold J. E., Lambert-Lacroix S., 
#' Picard F. (2017). High Dimensional Classification with combined Adaptive 
#' Sparse PLS and Logistic Regression, (in prep), 
#' available on (\url{http://arxiv.org/abs/1502.05933}).
#' 
#' Chun, H., & Keles, S. (2010). Sparse partial least squares regression for 
#' simultaneous dimension reduction and variable selection.  Journal of the 
#' Royal Statistical Society. Series B (Methodological), 72(1), 3-25. 
#' doi:10.1111/j.1467-9868.2009.00723.x
#' 
#' Chung, D., & Keles, S. (2010). Sparse partial least squares classification 
#' for high dimensional data. Statistical Applications in Genetics and 
#' Molecular Biology, 9, Article17. doi:10.2202/1544-6115.1492
#' 
#' @author
#' Ghislain Durif (\url{http://thoth.inrialpes.fr/people/gdurif/}). 
#' 
#' Adapted in part from spls code by H. Chun, D. Chung and S.Keles 
#' (\url{http://cran.r-project.org/web/packages/spls/index.html}).
#' 
#' @seealso \code{\link{logit.spls.cv}}
#' 
#' @examples
#' \dontrun{
#' ### load plsgenomics library
#' library(plsgenomics)
#' 
#' ### generating data
#' n <- 100
#' p <- 1000
#' sample1 <- sample.cont(n=100, p=1000, kstar=20, lstar=2, 
#'                        beta.min=0.25, beta.max=0.75, mean.H=0.2, 
#'                        sigma.H=10, sigma.F=5, sigma.E=5)
#'                        
#' X <- sample1$X
#' Y <- sample1$Y
#' 
#' ### tuning the hyper-parameters
#' cv1 <- spls.adapt.tune(X=X, Y=Y, lambda.l1.range=seq(0.05, 0.95, by=0.1), 
#'                        ncomp.range=1:2, weight.mat=NULL, adapt=TRUE, 
#'                        center.X=TRUE, center.Y=TRUE, 
#'                        scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, 
#'                        return.grid=TRUE, ncores=1, nfolds=10, nrun=1)
#' str(cv1)
#' 
#' ### otpimal values
#' cv1$lambda.l1.opt
#' cv1$ncomp.opt
#' }
#' 
#' @export
spls.cv <- function(X, Y, lambda.l1.range, ncomp.range, weight.mat=NULL, 
                    adapt=TRUE, center.X=TRUE, center.Y=TRUE, 
                    scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, 
                    return.grid=FALSE, ncores=1, nfolds=10, nrun=1) {
	
	#####################################################################
	#### Initialisation
	#####################################################################
	X <- as.matrix(X)
	n <- nrow(X) # nb observations
	p <- ncol(X) # nb covariates
	index.p <- c(1:p)
	Y <- as.matrix(Y)
	q = ncol(Y)
	one <- matrix(1,nrow=1,ncol=n)
	
	#####################################################################
	#### Tests on type input
	#####################################################################
	
	# On weighting matrix V
	if(!is.null(weight.mat)) { # weighting in scalar product (in observation space of dimension n)
		Vfull <- as.matrix(weight.mat) 
		
		if ((!is.matrix(Vfull)) || (!is.numeric(Vfull))) {
			stop("Message from spls.adapt: Vfull is not of valid type")}
		
		if ((n != ncol(Vfull)) || (n != nrow(Vfull))) {
			stop("Message from spls.adapt: wrong dimension for Vfull, must be a square matrix of size the number of observations in Xtrain")
		}
	} else { # no weighting in sclar product
		Vfull <- diag(rep(1,n), nrow=n, ncol=n)
	}
	
	# On weighted.center
	if ( (weighted.center) && (is.null(weight.mat))) {
		stop("Message from spls.adapt: if the centering is weighted, the weighting matrix V should be provided")
	}
	
	# ncores
	if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
		stop("Message from spls.cv: ncores is not of valid type")
	}
	
	# nfolds
	if ((!is.numeric(nfolds)) || (round(nfolds)-nfolds!=0) || (nfolds<1)) {
		stop("Message from spls.cv: nfolds is not of valid type")
	}
	
	
	#####################################################################
	#### Cross-validation: computation on each fold over the entire grid
	#####################################################################
	
	## the train set is partitioned into nfolds part, each observation is assigned into a fold
	fold.obs <- sort(rep(1:nfolds, length.out = n))
	
	## hyper-parameter grid
	grid <- expand.grid(lambda.l1=lambda.l1.range, ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
	
	cv.grid.allfolds <- matrix( unlist( mclapply(1:nfolds, function(k) {
		
		
		#### train and test variable
		Xtrain <- subset(X, fold.obs != k)
		Ytrain <- subset(Y, fold.obs != k)
		
		ntrain <- nrow(Xtrain)
		
		Xtest <- subset(X, fold.obs == k)
		Ytest <- subset(Y, fold.obs == k)
		
		ntest <- nrow(Xtest)
		
		V <- Vfull[fold.obs == k, fold.obs == k]
		
		#####################################################################
		#### Tests on type input
		#####################################################################
		# On Xtrain
		if ((!is.matrix(Xtrain)) || (!is.numeric(Xtrain))) {
			stop("Message from spls.adapt: Xtrain is not of valid type")
		}
		
		if (p==1) {
			stop("Message from spls.adapt: p=1 is not valid")}
		
		# On Xtest if necessary	
		if (is.vector(Xtest)==TRUE) {
			Xtest <- matrix(Xtest,nrow=1)
		}
		
		Xtest <- as.matrix(Xtest)
		ntest <- nrow(Xtest) 
		
		if ((!is.matrix(Xtest)) || (!is.numeric(Xtest))) {
			stop("Message from spls.adapt: Xtest is not of valid type")}
		
		if (p != ncol(Xtest)) {
			stop("Message from spls.adapt: columns of Xtest and columns of Xtrain must be equal")
		}
		
		# On Ytrain
		if ((!is.matrix(Ytrain)) || (!is.numeric(Ytrain))) {
			stop("Message from spls.adapt: Ytrain is not of valid type")
		}
		
		if (q != 1) {
			stop("Message from spls.adapt: Ytrain must be univariate")
		}
		
		if (nrow(Ytrain)!=ntrain) {
			stop("Message from spls.adapt: the number of observations in Ytrain is not equal to the Xtrain row number")
		}
		
		# On hyper parameter: lambda.ridge, lambda.l1
		if ( (sum(!is.numeric(lambda.l1.range))) || (sum(lambda.l1.range<0)) || (sum(lambda.l1.range>1)) ) {
			stop("Message from spls.adapt: lambda is not of valid type")
		}
		
		# ncomp type
		if ( (sum(!is.numeric(ncomp.range))) || (sum(round(ncomp.range)-ncomp.range!=0)) || (sum(ncomp.range<1)) || (sum(ncomp.range>p)) ) {
			stop("Message from spls.adapt: ncomp is not of valid type")
		}
		
		
		#####################################################################
		#### centering and scaling
		#####################################################################
		if (!weighted.center) {
			
			# Xtrain mean
			meanXtrain <- apply(Xtrain, 2, mean)
			
			# Xtrain sd
			sigmaXtrain <- apply(Xtrain, 2, sd)
			# test if predictors with null variance
			if ( any( sigmaXtrain < .Machine$double.eps )) {
				stop("Some of the columns of the predictor matrix have zero variance.")
			}
			
			# centering & eventually scaling X
			if(center.X && scale.X) {
				sXtrain <- scale( Xtrain, center=meanXtrain, scale=sigmaXtrain)
			} else if(center.X && !scale.X) {
				sXtrain <- scale( Xtrain, center=meanXtrain, scale=FALSE)
			} else {
				sXtrain <- Xtrain
			}
			
			# Y mean
			meanYtrain <- apply(Ytrain, 2, mean)
			
			# Y sd
			sigmaYtrain <- apply(Ytrain, 2, sd)
			# test if predictors with null variance
			if ( any( sigmaYtrain < .Machine$double.eps )) {
				stop("The response matrix has zero variance.")
			}
			# centering & eventually scaling Y
			if(center.Y && scale.Y) {
				sYtrain <- scale( Ytrain, center=meanYtrain, scale=sigmaYtrain)
			} else if(center.Y && !scale.Y) {
				sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE)
			} else {
				sYtrain <- Ytrain
			}
			
			if(center.Y && scale.Y) {
				sYtest <- scale( Ytest, center=meanYtrain, scale=sigmaYtrain)
			} else if(center.Y && !scale.Y) {
				sYtest <- scale( Ytest, center=meanYtrain, scale=FALSE)
			} else {
				sYtest <- Ytest
			}
			
			# Xtest	
			## centering and scaling depend on Xtest
			if(center.X && scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=sigmaXtrain )
			} else if(center.X && !scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
			} else {
				sXtest <- Xtest
			}
			
		} else { # weighted scaling
			
			sumV <- sum(diag(V))
			
			# X mean
			meanXtrain <- matrix(diag(V), nrow=1) %*% Xtrain / sumV
			
			# X sd
			sigmaXtrain <- apply(Xtrain, 2, sd)
			# test if predictors with null variance
			if ( any( sigmaXtrain < .Machine$double.eps ) ) {
				stop("Some of the columns of the predictor matrix have zero variance.")
			}
			# centering & eventually scaling X
			sXtrain <- scale( Xtrain, center=meanXtrain, scale=FALSE )
			
			# Y mean
			meanYtrain <- matrix(diag(V), nrow=1) %*% Ytrain / sumV
			
			# Y sd
			sigmaYtrain <- apply(Ytrain, 2, sd)
			# test if predictors with null variance
			if ( any( sigmaYtrain < .Machine$double.eps ) ) {
				stop("The response matrix have zero variance.")
			}
			
			# centering & eventually scaling Y
			sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE )
			sYtest <- scale( Ytest, center=meanYtrain, scale=FALSE )
			
			# Xtest
			sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
			
		}
		
		
		#####################################################################
		#### Computation over the entire grid
		#####################################################################
		
		### sur chaque fold, on calcule pour tout lambda.ridge.range
		cv.grid.byfold <- matrix( sapply( split(grid, f=row.names(grid)), function(grid.line) {
			
			model <- tryCatch( spls.aux(Xtrain=Xtrain, sXtrain=sXtrain, 
			                            Ytrain=Ytrain, sYtrain=sYtrain, 
			                            lambda.l1=grid.line$lambda.l1, 
			                            ncomp=grid.line$ncomp, 
			                            weight.mat=weight.mat, 
			                            Xtest=Xtest, sXtest=sXtest, 
			                            adapt=adapt, meanXtrain=meanXtrain, 
			                            meanYtrain=meanYtrain, 
			                            sigmaXtrain=sigmaXtrain, 
			                            sigmaYtrain=sigmaYtrain, 
			                            center.X=center.X, center.Y=center.Y,
			                            scale.X=scale.X, scale.Y=scale.Y, 
			                            weighted.center=weighted.center), 
			                   error = function(e) { warnings("Message from spls.adapt.cv: error when fitting a model in crossvalidation"); return(NULL);} )
			
			## resutls
			res <- numeric(4)
			
			if(!is.null(model)) {
				res <- c(grid.line$lambda.l1, grid.line$ncomp, k, sum((model$hatYtest - sYtest)^2) / ntest)
			} else {
				res <- c(grid.line$lambda.l1, grid.line$ncomp, k, NA)
			}
			
			return(res)
			
		}), ncol=4, byrow=TRUE)
		
		return( t(cv.grid.byfold) )
		
		
	}, mc.cores = ncores, mc.silent=TRUE)), ncol=4, byrow=TRUE)
	
	cv.grid.allfolds <- data.frame(cv.grid.allfolds)
	colnames(cv.grid.allfolds) <- c("lambda.l1", "ncomp", "nfold", "error")
	
	#####################################################################
	#### Find the optimal point in the grid
	#####################################################################
	
	
	## compute the mean error over the folds for each point of the grid
	cv.grid <- data.frame( as.matrix( with( cv.grid.allfolds, aggregate(error, list(lambda.l1, ncomp), function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)) ))))
	colnames(cv.grid) <- c("lambda.l1", "ncomp", "error", "error.sd")
	
	
	##### return
	if(return.grid) {
		
		return( list(lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], cv.grid=cv.grid) )
		
	} else {
		
		return( list(lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], cv.grid=NULL) )
		
	}
	
	
}
