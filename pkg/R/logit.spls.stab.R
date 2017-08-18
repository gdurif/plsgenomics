### logit.spls.stab.R  (2014-10)
###
###    Tuning parameters (ncomp, lambda.l1, lambda.ridge) for Ridge Iteratively Reweighted Least Squares followed by Adaptive Sparse PLS regression for binary response, by K-fold cross-validation
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


logit.spls.stab <- function(X, Y, lambda.ridge.range, lambda.l1.range, 
                            ncomp.range, 
                            adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                            ncores=1, nresamp=100, 
                            center.X=TRUE, scale.X=FALSE, weighted.center=TRUE, 
                            seed=NULL, verbose=TRUE) {
     
     #####################################################################
     #### Initialisation
     #####################################################################
     X <- as.matrix(X)
     n <- nrow(X) # nb observations
     p <- ncol(X) # nb covariates
     index.p <- c(1:p)
     if(is.factor(Y)) {
          Y <- as.numeric(levels(Y))[Y]
     }
     Y <- as.integer(Y)
     Y <- as.matrix(Y)
     q <- ncol(Y)
     one <- matrix(1,nrow=1,ncol=n)
     
     if(!is.null(seed)) {
          set.seed(seed)
     }
     
     cnames <- NULL
     if(!is.null(colnames(X))) {
          cnames <- colnames(X)
     } else {
          cnames <- paste0(1:p)
     }
     
     
     #####################################################################
     #### Tests on type input
     #####################################################################
     
     # if multicategorical response
     if(length(table(Y)) > 2) {
          warning("message from logit.spls.stab: multicategorical response")
          results = m.logit.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, lambda.l1.range=lambda.l1.range, 
                                      ncomp.range=ncomp.range, 
                                      adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
                                      ncores=ncores, nresamp=nresamp, nfolds=nfolds, 
                                      center.X=center.X, scale.X=scale.X, weighted.center=weighted.center, 
                                      seed=seed, verbose=verbose)
          return(results)
     }
     
     # On X
     if ((!is.matrix(X)) || (!is.numeric(X))) {
          stop("message from logit.spls.stab: X is not of valid type")
     }
     
     if (p==1) {
          stop("message from logit.spls.stab: p=1 is not valid")
     }
     
     # On Y
     if ((!is.matrix(Y)) || (!is.numeric(Y))) {
          stop("message from logit.spls.stab: Y is not of valid type")
     }
     
     if (q != 1) {
          stop("message from logit.spls.stab: Y must be univariate")
     }
     
     if (nrow(Y)!=n) {
          stop("message from logit.spls.stab: the number of observations in Y is not equal to the number of row in X")
     }
     
     # On Y value
     if (sum(is.na(Y))!=0) {
          stop("message from logit.spls.stab: NA values in Ytrain")
     }
     
     if (sum(!(Y %in% c(0,1)))!=0) {
          stop("message from logit.spls.stab: Y is not of valid type")
     }
     
     if (sum(as.numeric(table(Y))==0)!=0) {
          stop("message from logit.spls.stab: there are empty classes")
     }
     
     # On hyper parameter: lambda.ridge, lambda.l1
     if (any(!is.numeric(lambda.ridge.range)) || any(lambda.ridge.range<0) 
         || any(!is.numeric(lambda.l1.range)) || any(lambda.l1.range<0) || any(lambda.l1.range>1)) {
          stop("Message from logit.spls.stab: lambda is not of valid type")
     }
     
     # ncomp type
     if (any(!is.numeric(ncomp.range)) || any(round(ncomp.range)-ncomp.range!=0) || any(ncomp.range<0) || any(ncomp.range>p)) {
          stop("Message from logit.spls.stab: ncomp is not of valid type")
     }
     
     # maxIter
     if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
          stop("message from logit.spls.stab: maxIter is not of valid type")
     }
     
     # ncores
     if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
          stop("message from logit.spls.stab: ncores is not of valid type")
     }
     
     # necessary to insure that both classes are represented in each folds
     if( any(as.vector(table(Y))<nfolds)) {
          stop("message from logit.spls.stab: there is a class defined by Y that has less members than the number of folds nfold")
     }
     
     
     #####################################################################
     #### Stability selection procedure
     #####################################################################
     
     ## computation on the folds x run grid
     grid.resampling <- as.matrix( Reduce("rbind", mclapply(1:nresamp, function(id.samp) {
          
          #### train and test variable
          ntrain = floor(0.5*n)
          ntest = n - ntrain
          
          index.train = sort(sample(1:n, size=ntrain))
          index.test = (1:n)[-index.train]
          
          Xtrain = X[index.train,]
          Ytrain = Y[index.train]
          
          Xtest = X[index.test,]
          Ytest = Y[index.test]
          
          condition = any(table(Ytrain)<1)
          test = 0
          while(condition & test<100) {
               index.train = sort(sample(1:n, size=ntrain))
               index.test = (1:n)[-index.train]
               
               Xtrain = X[index.train,]
               Ytrain = Y[index.train]
               
               Xtest = X[index.test,]
               Ytest = Y[index.test]
               
               condition = any(table(Ytrain)<1)
               test = test+1
          }
          
          if(test==100) {
               stop("message from logit.spls.stab: error in subsampling, a sub-sample has only one Y-label.")
          }
          
          #### hyper-parameter grid
          paramGrid <- expand.grid(lambdaL1=lambda.l1.range, 
                                   lambdaL2=lambda.ridge.range, 
                                   ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
          
          #### fit the model for the different lambda.l1
          grid_out <- as.matrix( Reduce("rbind", lapply(split(paramGrid, f=row.names(paramGrid)), function(gridRow) {
               
               lambdaL1 <- gridRow$lambdaL1
               lambdaL2 <- gridRow$lambdaL2
               ncomp <- gridRow$ncomp
               
               ## fit the model for the chosen lambda ridge
               fit_out <- logit.spls(Xtrain=Xtrain, Ytrain=Ytrain, 
                                     lambda.ridge=lambdaL2, 
                                     lambda.l1=lambdaL1, 
                                     ncomp=ncomp, Xtest=Xtest, adapt=adapt, 
                                     maxIter=maxIter, svd.decompose=svd.decompose, 
                                     center.X=center.X, scale.X=scale.X, 
                                     weighted.center=weighted.center)
               
               ## selected variables ?
               sel_var <- fit_out$Anames
               status_var <- rep(0, length(cnames))
               status_var[which(cnames %in% sel_var)] <- rep(1, length(which(cnames %in% sel_var)))
               
               tmp <- c(lambdaL1, lambdaL2, ncomp, id.samp, sum(status_var), status_var)
               return(tmp)
               
          })))
          
          rownames(grid_out) <- NULL
          return(grid_out)
          
     }, mc.cores=ncores, mc.silent=!verbose)))
     
     grid.resampling <- data.frame(grid.resampling)
     colnames(grid.resampling) <- c("lambdaL1", "lambdaL2", "ncomp", 
                                    "id", "nbVar", cnames)
     
     grid.resampling$point = paste0(grid.resampling$lambdaL1, "_",
                                    grid.resampling$lambdaL2, "_",
                                    grid.resampling$ncomp)
     
     o.grid <- order(grid.resampling$nbVar)
     grid.resampling <- grid.resampling[o.grid,]
     
     if(any(table(grid.resampling$point)<nresamp)) {
          warning("message from logit.spls.stab: empty classe in a resampling")
          print(table(grid.resampling$point))
     }
     
     #####################################################################
     #### Compute q_lambda
     #####################################################################
     
     ## increasing value of q_lambda
     tmp_qLambda <- as.matrix( Reduce("rbind", mclapply(1:nresamp, function(id.samp) {
          
          tmp1 <- subset(grid.resampling, id==id.samp)
          tmp2 <- apply(t(tmp1)[-c(1:5,tail(1:ncol(grid.resampling), 1)),],1,cumsum) # cumsum by genes
          tmp3 <- apply(t(tmp2), 2, function(x) return(sum(x!=0)))
          
          tmp4 <- cbind(tmp1[,c(1:4)], unname(tmp3))
          
          return(tmp4)
          
     }, mc.cores=ncores, mc.silent=!verbose)))
     tmp_qLambda <- data.frame(tmp_qLambda)
     colnames(tmp_qLambda) <- c("lambdaL1", "lambdaL2", "ncomp", "id.samp", "qLambda")
     
     qLambda <- ddply(tmp_qLambda, c("lambdaL1", "lambdaL2", "ncomp"), 
                      function(x) colMeans(x[c("qLambda")], na.rm=TRUE))
     o.qLambda <- order(qLambda$qLambda)
     qLambda <- qLambda[o.qLambda,]
     
     #####################################################################
     #### Compute p_j(lambda)
     #####################################################################
     
     probs_lambda <- ddply(grid.resampling, c("lambdaL1", "lambdaL2", "ncomp"), 
                           function(x) colMeans(x[cnames], na.rm=TRUE))
     probs_lambda <- probs_lambda[o.qLambda,]
     
     # select variables in another function
     
     
     ##### return
     return(list(qLambda=qLambda, probs.lambda=probs_lambda))
     
}


