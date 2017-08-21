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


stability.selection <- function(stab.out, piThreshold=0.6, rhoError=10) {
     
     qLambda <- stab.out$q.Lambda
     probs_lambda <- stab.out$probs.lambda
     
     p <- ncol(probs_lambda[,-1])
     
     ## select variables
     qLambdaMax = sqrt((2*piThreshold-1) * p * rhoError)
     
     lambda_ok <- qLambda$lambda[which(qLambda$qLambda <= qLambdaMax)]
     
     tmp_probs <- probs_lambda[probs_lambda$lambda %in% lambda_ok,]
     
     which_var <-  apply(tmp_probs[,-1], 2, function(x) max(x) >= piThreshold)
     
     selected.variables <- colnames(probs_lambda[,-1])[which_var]
     
     return(selected.variables)
}


stability.selection.heatmap <- function(stab_out, ...) {
     matrix.heatmap(as.matrix(stab_out$probs_lambda[,-c(1:3)]),
                    xlab="covariates", ylab="hyper-param. values",
                    ...)
}
