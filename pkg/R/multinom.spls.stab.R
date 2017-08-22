### rirls.spls.stab.R  (2014-10)
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


m.rirls.spls.stab <- function(X, Y, lambda.ridge.range, lambda.l1.range, ncomp=1, 
                            adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                            ncores=1, nresamp=100, nfolds=5, 
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
     
     # if Binary response
     if(length(table(Y)) == 2) {
          warning("message from m.rirls.spls.tune: binary response")
          results = rirls.spls.stab(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, lambda.l1.range=lambda.l1.range, ncomp=ncomp, 
                                    adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
                                    ncores=ncores, nresamp=nresamp, nfolds=nfolds, 
                                    center.X=center.X, scale.X=scale.X, weighted.center=weighted.center, 
                                    seed=seed, verbose=verbose)
          return(results)
     }
     
     # On X
     if ((!is.matrix(X)) || (!is.numeric(X))) {
          stop("Message from m.rirls.spls.tune: X is not of valid type")
     }
     
     if (p==1) {
          # stop("Message from m.rirls.spls.tune: p=1 is not valid")
          warning("Message from m.rirls.spls.tune: p=1 is not valid, ncomp.range is set to 0")
          ncomp.range <- 0
     }
     
     # On Y
     if ((!is.matrix(Y)) || (!is.numeric(Y))) {
          stop("Message from m.rirls.spls: Y is not of valid type")
     }
     
     if (q != 1) {
          stop("Message from m.rirls.spls: Y must be univariate")
     }
     
     if (nrow(Y)!=n) {
          stop("Message from m.rirls.spls: the number of observations in Y is not equal to the Xtrain row number")
     }
     
     # On Y value
     if (sum(is.na(Y))!=0) {
          stop("Message from m.rirls.spls: NA values in Y")
     }
     
     if((sum(floor(Y)-Y)!=0)||(sum(Y<0)>0)) {
          stop("Message from m.rirls.spls: Y is not of valid type")
     }
     
     if (sum(as.numeric(table(Y))==0)!=0) {
          stop("Message from m.rirls.spls: there are empty classes")
     }
     
     
     # maxIter
     if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
          stop("message from rirls.spls.stab: maxIter is not of valid type")
     }
     
     # ncores
     if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
          stop("message from rirls.spls.stab: ncores is not of valid type")
     }
     
     # nfolds
     if ((!is.numeric(nfolds)) || (round(nfolds)-nfolds!=0) || (nfolds<1)) {
          stop("message from rirls.spls.stab: nfolds is not of valid type")
     }
     # necessary to insure that both classes are represented in each folds
     if( any(as.vector(table(Y))<nfolds)) {
          stop("message from rirls.spls.stab: there is a class defined by Y that has less members than the number of folds nfold")
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
          
          condition = any(table(Ytrain)<nfolds)
          test = 0
          while(condition & test<100) {
               index.train = sort(sample(1:n, size=ntrain))
               index.test = (1:n)[-index.train]
               
               Xtrain = X[index.train,]
               Ytrain = Y[index.train]
               
               Xtest = X[index.test,]
               Ytest = Y[index.test]
               
               condition = any(table(Ytrain)<nfolds)
               test = test+1
          }
          
          if(test==100) {
               return(NULL)
          }
          
          #### fit the model for the different lambda.l1
          grid_out <- sapply(lambda.l1.range, function(lambda) {
               
               ## tune the lambda ridge param
               cv_out <- m.rirls.spls.tune(X=Xtrain, Y=Ytrain, lambda.ridge.range=lambda.ridge.range, 
                                           lambda.l1.range=lambda, ncomp.range=ncomp, adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
                                           return.grid=FALSE, ncores=1, nfolds=nfolds, nrun=1, 
                                           center.X=center.X, scale.X=scale.X, weighted.center=weighted.center, 
                                           seed=NULL, verbose=FALSE)
               
               ## fit the model for the chosen lambda ridge
               fit_out <- m.rirls.spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.ridge=cv_out$lambda.ridge.opt, 
                                       lambda.l1=lambda, 
                                       ncomp=ncomp, Xtest=Xtest, adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
                                       center.X=center.X, scale.X=scale.X, weighted.center=weighted.center)
               
               ## selected variables ?
               sel_var <- fit_out$Anames
               status_var <- rep(0, length(cnames))
               status_var[which(cnames %in% sel_var)] <- rep(1, length(which(cnames %in% sel_var)))
               
               tmp <- c(lambda, id.samp, status_var)
               return(tmp)
               
          })
          
          return(t(grid_out))
          
     }, mc.cores=ncores, mc.silent=!verbose)))
     
     grid.resampling <- data.frame(grid.resampling)
     colnames(grid.resampling) <- c("lambda", "id", cnames)
     
     o.lambda <- order(grid.resampling$lambda, decreasing=TRUE)
     grid.resampling <- grid.resampling[o.lambda,]
     
     if(any(table(grid.resampling$lambda)<nresamp)) {
          warning("message from rirls.spls.stab: empty classe in a resampling")
          print(table(grid.resampling$lambda))
     }
     
     #####################################################################
     #### Compute q_lambda
     #####################################################################
     
     tmp_qLambda <- as.matrix( Reduce("rbind", mclapply(1:nresamp, function(id.samp) {
          
          tmp1 <- subset(grid.resampling, id==id.samp)
          tmp2 <- apply(t(tmp1)[-c(1,2),],1,cumsum) # cumsum by genes
          tmp3 <- apply(t(tmp2), 2, function(x) return(sum(x!=0)))
          
          tmp4 <- cbind(sort(lambda.l1.range, decreasing=TRUE), rep(id.samp, length(lambda.l1.range)), unname(tmp3))
          
          return(tmp4)
          
     }, mc.cores=ncores, mc.silent=!verbose)))
     tmp_qLambda <- data.frame(tmp_qLambda)
     colnames(tmp_qLambda) <- c("lambda", "id", "qLambda")
     
     qLambda <- ddply(tmp_qLambda, "lambda", function(x) colMeans(x[c("qLambda")], na.rm=TRUE))
     o.lambda2 <- order(qLambda$lambda, decreasing=TRUE)
     qLambda <- qLambda[o.lambda2,]
     
     #####################################################################
     #### Compute p_j(lambda)
     #####################################################################
     
     probs_lambda <- ddply(grid.resampling, "lambda", function(x) colMeans(x[cnames], na.rm=TRUE))
     probs_lambda <- probs_lambda[o.lambda2,]
     
     # select variables in another function
     
     
     ##### return
     return(list(qLambda=qLambda, probs_lambda=probs_lambda))
     
}


