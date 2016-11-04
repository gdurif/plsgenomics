# load plsgenomics library
library(plsgenomics)

# load the Ecoli data
data(Ecoli)

# perform pls regression
# with unit latent components
X = Ecoli$CONNECdata
Y = Ecoli$GEdata
K = 3
pls.regression(Xtrain=X,Ytrain=Y,Xtest=Ecoli$CONNECdata,
               ncomp=1:K,unit.weights=FALSE)


var.rec = function(model, X, Y, K) {
     
     # data
     X = as.matrix(X)
     Y = as.matrix(Y)
     
     X = scale(X, center=TRUE, scale=FALSE)
     Y = scale(Y, center=TRUE, scale=FALSE)
     
     n = nrow(X)
     p = ncol(X)
     q = ncol(Y)
     
     T.K = matrix(model$T, nrow=n, ncol=K)
     
     # stockage projections sur T_1,..., T_k
     X.res = matrix(NA, nrow = n, ncol = p)
     Y.res = matrix(NA, nrow = n, ncol = q)
     
     res <- sapply(1:K, function(k) {
          
          T.k <- as.matrix(T.K[,1:k])
          
          # projection of X onto T_1, ..., T_k
          X.res = T.k %*% solve(t(T.k) %*% T.k) %*% t(T.k) %*% X
          var.X1 = sum(diag( var(X.res)))/sum(diag( var(X)))
          
          # projection of Y onto T_1, ..., T_k
          Y.res = T.k %*% solve(t(T.k) %*% T.k) %*% t(T.k) %*% Y
          var.Y1 = sum(diag( var(Y.res)))/sum(diag( var(Y)))
          
          T.k <- as.matrix(T.K[,k])
          
          # projection of X onto T_k
          X.res = T.k %*% solve(t(T.k) %*% T.k) %*% t(T.k) %*% X
          var.X2 = sum(diag( var(X.res)))/sum(diag( var(X)))
          
          # projection of Y onto T_k
          Y.res = T.k %*% solve(t(T.k) %*% T.k) %*% t(T.k) %*% Y
          var.Y2 = sum(diag( var(Y.res)))/sum(diag( var(Y)))
          
          return(c(var.X1, var.X2, var.Y1, var.Y2))
     })
     
     res <- data.frame(t(res))
     colnames(res) = c("cum.exp.var.X", "exp.var.X", "cum.exp.var.Y", "exp.var.Y")
     rownames(res) = paste0("comp", 1:K)
     
     return(res)
     
}

res = var.rec(model, Ecoli$CONNECdata, Ecoli$GEdata, 3)
res
cbind(res[,1], cumsum(res[,2]))

vip = function(model, X, Y, K) {
     
     # data
     X = as.matrix(X)
     Y = as.matrix(Y)
     
     X = scale(X, center=TRUE, scale=FALSE)
     Y = scale(Y, center=TRUE, scale=FALSE)
     
     n = nrow(X)
     p = ncol(X)
     q = ncol(Y)
     
     T.K = matrix(model$T, nrow=n, ncol=K)
     W.K = matrix(model$R, nrow=p, ncol=K)
     
     ### VIP 
     VIP = matrix(NA, nrow=p, ncol=K)
     cor2 = cor(Y, T.K, use="pairwise")^2
     
     VIP[,1] = W.K[,1]^2
     if(K>1) {
          for(k in 2:K) {
               if(q == 1) {
                    Rd = as.matrix(cor2[,1:k])
                    VIP[,k] = ( W.K[,1:k]^2 %*% Rd ) / sum(Rd)
               } else {
                    Rd = as.matrix(apply(cor2[,1:k], 2, sum))
                    VIP[,k] = ( W.K[,1:k]^2 %*% Rd ) / sum(Rd)
               }
          }
     }
     
     VIP = sqrt(p * VIP)
     if(!is.null(rownames(X))) {
          rownames(VIP) = colnames(X)
     } else {
          rownames(VIP) = as.character(1:p)
     }
     colnames(VIP) = paste0("comp", 1:K)
     
     return(VIP)
}

vip(model, Ecoli$CONNECdata, Ecoli$GEdata, 3)
