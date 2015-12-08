### rirls.spls.tune.R  (2015-10)
###
###    Tuning parameters (ncomp, lambda.l1, lambda.ridge) for Ridge Iteratively Reweighted 
###    Least Squares followed by Adaptive Sparse PLS regression for multicategorial response, 
###    by K-fold cross-validation
###
### Copyright 2015-10 Ghislain DURIF
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


m.rirls.spls.tune <- function(X, Y, lambda.ridge.range, lambda.l1.range, ncomp.range, adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                            return.grid=FALSE, ncores=1, nfolds=10, nrun=1, 
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
     
     
     #####################################################################
     #### Tests on type input
     #####################################################################
     
     # On X
     if ((!is.matrix(X)) || (!is.numeric(X))) {
          stop("Message from m.rirls.spls.tune: X is not of valid type")
     }
     
     if (p==1) {
          stop("Message from m.rirls.spls.tune: p=1 is not valid")
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
     
     # On hyper parameter: lambda.ridge, lambda.l1
     if ( any(!is.numeric(lambda.ridge.range)) || any(lambda.ridge.range<0) || any(!is.numeric(lambda.l1.range)) || any(lambda.l1.range<0)) {
          stop("Message from m.rirls.spls: lambda is not of valid type")
     }
     
     # ncomp type
     if ( any(!is.numeric(ncomp.range)) || any(round(ncomp.range)-ncomp.range!=0) || any(ncomp.range<0) || any(ncomp.range>p)) {
          stop("Message from m.rirls.spls: ncomp.range is not of valid type")
     }
     
     
     # maxIter
     if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
          stop("Message from m.rirls.spls.tune: maxIter is not of valid type")
     }
     
     # ncores
     if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
          stop("Message from m.rirls.spls.tune: ncores is not of valid type")
     }
     
     # nfolds
     if ((!is.numeric(nfolds)) || (round(nfolds)-nfolds!=0) || (nfolds<1)) {
          stop("Message from m.rirls.spls.tune: nfolds is not of valid type")
     }
     
     # necessary to insure that both classes are represented in each folds
     if( any(as.vector(table(Y))<nfolds) ) {
          stop("Message from m.rirls.spls.tune: there is a class defined by Y that has less members than the number of folds nfold")
     }
     
     ## fold size
     fold.size = n %/% nfolds
     if( length(unique(Y))>fold.size ) {
          stop("Message from m.rirls.spls.tune: there are too many classes compared to the size of folds")
     }
     
     
     #####################################################################
     #### Cross-validation: computation on each fold over the entire grid
     #####################################################################
     
     ## draw nrun nfolds partition
     folds.obs = sapply(1:nrun, function(run) {
          ## the train set is partitioned into nfolds part, each observation is assigned into a fold
          
          # observations from the different classes
          classes <- sort(unique(Y))
          index <- lapply(classes, function(g) {
               return((1:n)[Y==g])
          })
          
          # each folds contains at least one elements of each class
          # choose the represent of each class in each fold
          folds <- lapply(classes+1, function(g) {
               return(c(1:nfolds, rep(0, length(index[[g]]) - nfolds))[sample(x=1:length(index[[g]]), size=length(index[[g]]), replace=FALSE)])
          })
          
          ind.folds <- lapply(classes+1, function(g) {
               return((index[[g]])[folds[[g]]!=0])
          })
          
          # group the remaining observation and equally dispatch them into the folds
          rest <- unlist(lapply(classes+1, function(g) {
               return((index[[g]])[folds[[g]]==0])
          }))
          
          folds.rest = rep(1:nfolds, length.out = length(rest))[sample(x=1:length(rest), size=length(rest), replace=FALSE)]
          
          folds.res = rep(NA, n)
          
          for(g in (classes+1)) {
               folds.res[ind.folds[[g]]] <- (folds[[g]])[folds[[g]]!=0]
          }
          folds.res[rest] = folds.rest
          
          #table(folds.res)
          #sapply(1:nfolds, function(x) table(Y[folds.res==x]))
          
          return(folds.res)
     })
     #for(i in 1:nfolds) print(table(Y[folds.obs[,1]==i,]))
     
     
     
     ## hyper-parameter grid
     grid = expand.grid(lambda.ridge=lambda.ridge.range, lambda.l1=lambda.l1.range, ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
     
     ## folds x run grid
     folds.grid = expand.grid(k=1:nfolds, run=1:nrun, KEEP.OUT.ATTRS=FALSE)
     
     ## computation on the folds x run grid
     cv.grid.allfolds = matrix( unlist( mclapply(split(folds.grid, f=row.names(folds.grid)), function(folds.line) {
          
          
          #### train and test variable
          Xtrain <- subset(X, folds.obs[,folds.line$run] != folds.line$k)
          Ytrain <- subset(Y, folds.obs[,folds.line$run] != folds.line$k)
          
          ntrain <- nrow(Xtrain)
          
          Xtest <- subset(X, folds.obs[,folds.line$run] == folds.line$k)
          Ytest <- subset(Y, folds.obs[,folds.line$run] == folds.line$k)
          
          ntest <- nrow(Xtest)
          
          
          #####################################################################
          #### Move into the reduced space
          #####################################################################
          
          r <- p #min(p, ntrain)
          G <- max(Ytrain)
          DeletedCol <- NULL
          
          ### Standardize the Xtrain matrix
          # standard deviation (biased one) of Xtrain
          sigma2train <- apply(Xtrain, 2, var) * (ntrain-1)/(ntrain)
          
          # test on sigma2train
          # predictor with null variance ?
          if (sum(sigma2train < .Machine$double.eps)!=0){
               
               # predicteur with non null variance < 2 ?
               if (sum(sigma2train < .Machine$double.eps)>(p-2)){
                    stop("Message from m.rirls.spls: the procedure stops because number of predictor variables with no null variance is less than 1.")
               }
               
               warning("There are covariables with nul variance")
               
               # remove predictor with null variance
               Xtrain <- Xtrain[,which(sigma2train >= .Machine$double.eps)]
               if (!is.null(Xtest)) {
                    Xtest <- Xtest[,which(sigma2train>= .Machine$double.eps)]
               }
               
               # list of removed predictors
               DeletedCol <- index.p[which(sigma2train < .Machine$double.eps)]
               
               # removed null standard deviation
               sigma2train <- sigma2train[which(sigma2train >= .Machine$double.eps)]
               
               # new number of predictor
               p <- ncol(Xtrain)
               r <- p
          }
          
          # mean of Xtrain
          meanXtrain <- apply(Xtrain,2,mean)
          
          # center and scale Xtrain
          if(center.X && scale.X) {
               sXtrain <- scale(Xtrain, center=meanXtrain, scale=sqrt(sigma2train))
          } else if(center.X && !scale.X) {
               sXtrain <- scale(Xtrain, center=meanXtrain, scale=FALSE)
          } else {
               sXtrain <- Xtrain
          }
          
          sXtrain.nosvd <- sXtrain # keep in memory if svd decomposition
          
          # Compute the svd when necessary -> case p > ntrain (high dim)
          if ((p > ntrain) && (svd.decompose)) {
               # svd de sXtrain
               svd.sXtrain <- svd(t(sXtrain))
               # number of singular value non null
               r <- length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
               V <- svd.sXtrain$u[,1:r]
               D <- diag(c(svd.sXtrain$d[1:r]))
               U <- svd.sXtrain$v[,1:r]
               sXtrain <- U %*% D
               rm(D)
               rm(U)
               rm(svd.sXtrain)
          }
          
          # center and scale Xtest
          meanXtest <- apply(Xtest,2,mean)
          sigma2test <- apply(Xtest,2,var)
          
          if(center.X && scale.X) {
               sXtest <- scale(Xtest, center=meanXtrain, scale=sqrt(sigma2train))
          } else if(center.X && !scale.X) {
               sXtest <- scale(Xtest, center=meanXtrain, scale=FALSE)
          } else {
               sXtest <- Xtest
          }
          
          sXtest.nosvd <- sXtest # keep in memory if svd decomposition
          
          # if svd decomposition
          if ((p > ntrain) && (svd.decompose)) {
               sXtest <- sXtest%*%V
          }
          
          
          
          #####################################################################
          #### Computation over the entire grid
          #####################################################################
          
          ### sur chaque fold, on calcule pour tout lambda.ridge.range
          cv.grid.byfold <- matrix( sapply( split(grid, f=row.names(grid)), function(grid.line) {
               
               # On hyper parameter: lambda.ridge, lambda.l1
               if ((!is.numeric(grid.line$lambda.ridge)) || (grid.line$lambda.ridge<0) || (!is.numeric(grid.line$lambda.l1)) || (grid.line$lambda.l1<0)) {
                    stop("Message from m.rirls.spls.tune: lambda is not of valid type")
               }
               
               # ncomp type
               if ((!is.numeric(grid.line$ncomp)) || (round(grid.line$ncomp)-grid.line$ncomp!=0) || (grid.line$ncomp<0) || (grid.line$ncomp>p)) {
                    stop("Message from m.rirls.spls.tune: ncomp is not of valid type")
               }
               
               model <- tryCatch( m.rirls.spls.aux(sXtrain=sXtrain, sXtrain.nosvd=sXtrain.nosvd, Ytrain=Ytrain, lambda.ridge=grid.line$lambda.ridge, 
                                                 lambda.l1=grid.line$lambda.l1, ncomp=grid.line$ncomp, sXtest=sXtest, sXtest.nosvd=sXtest.nosvd, 
                                                 adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
                                                 meanXtrain=meanXtrain, sigma2train=sigma2train, 
                                                 center.X=center.X, scale.X=scale.X, weighted.center=weighted.center), 
                                  error = function(e) { print(e); warnings("Message from m.rirls.spls.tune: error when fitting a model in crossvalidation"); return(NULL);} )
               
               #model <- m.rirls.spls.aux(sXtrain=sXtrain, sXtrain.nosvd=sXtrain.nosvd, Ytrain=Ytrain, lambda.ridge=grid.line$lambda.ridge, 
               #                                    lambda.l1=grid.line$lambda.l1, ncomp=grid.line$ncomp, sXtest=sXtest, sXtest.nosvd=sXtest.nosvd, 
               #                                    adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, 
               #                                    meanXtrain=meanXtrain, sigma2train=sigma2train, 
               #                                    center.X=center.X, scale.X=scale.X, weighted.center=weighted.center)

               ## results
               res = numeric(6)
               
               if(!is.null(model)) {
                    res = c(grid.line$lambda.ridge, grid.line$lambda.l1, grid.line$ncomp, folds.line$k, folds.line$run, model$converged, sum(model$hatYtest != Ytest) / ntest)
               } else {
                    res = c(grid.line$lambda.ridge, grid.line$lambda.l1, grid.line$ncomp, folds.line$k, folds.line$run, NA, NA)
               }
               
               return(res)
               
          }), ncol=7, byrow=TRUE)
          
          return( as.vector(t(cv.grid.byfold)) )
          
          
     }, mc.cores=ncores, mc.silent=!verbose)), ncol=7, byrow=TRUE)
     #})), ncol=7, byrow=TRUE)
     
     cv.grid.allfolds = data.frame(cv.grid.allfolds)
     colnames(cv.grid.allfolds) = c("lambda.ridge", "lambda.l1", "ncomp", "nfold", "nrun", "converged", "error")
     
     #####################################################################
     #### Find the optimal point in the grid
     #####################################################################
     
     ## compute the number of NAs (=fails)
     cv.grid.fails <- data.frame( as.matrix( with( cv.grid.allfolds, aggregate(error, list(lambda.ridge, lambda.l1, ncomp), function(x) {sum(is.na(x))}))))
     colnames(cv.grid.fails) <- c("lambda.ridge", "lambda.l1", "ncomp", "nb.fail")
     
     ## check number of NAs (=fails)
     if(sum(cv.grid.fails$nb.fail>0.8*nrun*nfolds) > (0.8*nrow(grid))) {
          warnings("Message from m.rirls.spls.tune: too many errors during the cross-validation process, the grid is not enough filled")
     }     
     
     ## compute the mean error over the folds for each point of the grid
     cv.grid.error <- data.frame( as.matrix( with( cv.grid.allfolds, aggregate(error, list(lambda.ridge, lambda.l1, ncomp), function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)) ))))
     colnames(cv.grid.error) <- c("lambda.ridge", "lambda.l1", "ncomp", "error", "error.sd")
     
     ## compute the percentage of convergence over the folds for each point of the grid
     cv.grid.conv <- data.frame( as.matrix( with( cv.grid.allfolds, aggregate(converged, list(lambda.ridge, lambda.l1, ncomp), mean, na.rm=TRUE))))
     colnames(cv.grid.conv) <- c("lambda.ridge", "lambda.l1", "ncomp", "converged")
     
     ## % of convergence over all cross-validation process
     conv.per <- mean(cv.grid.conv$converged)
     
     ## merge the three tables
     cv.grid1 <- merge(cv.grid.error, cv.grid.conv, by = c("lambda.ridge", "lambda.l1", "ncomp"))
     cv.grid <- merge(cv.grid1, cv.grid.fails, by = c("lambda.ridge", "lambda.l1", "ncomp"))
     
     ##### return
     if(return.grid) {
          
          return( list(lambda.ridge.opt = cv.grid$lambda.ridge[which.min(cv.grid$error)], lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], conv.per=conv.per, cv.grid=cv.grid) )
          
     } else {
          
          return( list(lambda.ridge.opt = cv.grid$lambda.ridge[which.min(cv.grid$error)], lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], conv.per=conv.per, cv.grid=NULL) )
          
     }
     
  
}

