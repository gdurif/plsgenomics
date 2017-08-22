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


stability.selection.heatmap <- function(stab.out, ...) {
     matrix.heatmap(as.matrix(stab.out$probs.lambda[,-c(1:3)]),
                    xlab="covariates", ylab="selection sharpness",
                    ...)
}
