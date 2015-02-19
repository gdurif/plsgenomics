\name{spls.adapt}
\alias{spls.Adapt}

\title{Classification by Ridge Iteratively Reweighted Least Squares followed by Adaptive Sparse PLS regression for binary response}
\usage{
spls.adapt(Xtrain, Ytrain, lambda.l1, ncomp, weight.mat=NULL, Xtest=NULL, adapt=TRUE, 
           center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE)
}
\arguments{
	\item{Xtrain}{a (ntrain x p) data matrix of predictors. \code{Xtrain} must be a matrix. 
	Each row corresponds to an observation and each column to a predictor variable.}
	\item{Ytrain}{a ntrain vector of (continuous) responses. \code{Ytrain} must be a vector or a one column matrix. 
	and contains the response variable for each observation.}
	\item{lambda.l1}{a positive real value, in [0,1]. \code{lambda.l1} is the sparse penalty parameter 
	for the dimension reduction step by sparse PLS (see details).}
	\item{ncomp}{a positive integer. \code{ncomp} is the number of PLS components. 
	If \code{ncomp}=0,then the Ridge regression is performed without dimension reduction.}
	\item{weighted.mat}{a (ntrain x ntrain) matrix used to weight the l2 metric in observation spase if necessary,
	especially the covariance inverse of the Ytrain observations in heteroskedastic context. If NULL, the l2 metric
	is the standard one, corresponding to homoskedastic model.}
	\item{Xtest}{a (ntest x p) matrix containing the predictors for the test data set. \code{Xtest} may also be 
	a vector of length p (corresponding to only one test observation).If \code{Xtest} is not equal to NULL, then the 
	prediction step is made for these new predictor variables.}
	\item{adapt}{a boolean value, indicating whether the sparse PLS selection step sould be adaptive or nor.}
	\item{center.X}{a boolean value indicating whether the design matrices Xtrain in train set and Xtest in
	test set if non NULL should be centered or not}
	\item{scale.X}{a boolean value indicating whether the design matrices Xtrain in train set and Xtest in
	test set if non NULL should be scaled or not, scale.X=TRUE implies center.X=TRUE}
	\item{center.Y}{a boolean value indicating whether the response Ytrain in train set should be centered or not}
	\item{scale.Y}{a boolean value indicating whether the response Ytrain should be scaled or not, 
	scale.Y=TRUE implies center.Y=TRUE}
	\item{weighted.centered}{a boolean value indicating whether should the centering take into account the 
	weighted l2 metric or not (if TRUE, it implies that weighted.mat is non NULL).}
}
\description{
	The function \code{spls.adapt} performs compression, variable selection in regression context (with possible prediction) 
	using Durif et al. (2015) adaptive SPLS algorithm, based on sparse PLS developped by Chun and Kele\c{s} (2010).
  
}

\details{
	The columns of the data matrices \code{Xtrain} and \code{Xtest} may not be standardized, 
	since standardizing is can be performed by the function \code{spls.adapt} as a preliminary step
	before the algorithm is run. 
	
	The procedure described in Durif et al. (2015) is used to determine
	latent sparse components to be used for regression and when \code{Xtest} 
	is not equal to NULL, the procedure predicts the response for these new 
	predictor variables.  
 
}

Xtrain=Xtrain, Ytrain=Ytrain, sXtrain=sXtrain, sYtrain=sYtrain,
				 betahat=betahat, betahat.nc=betahat.nc,
				 meanXtrain=meanXtrain, meanYtrain=meanYtrain, sigmaXtrain=sigmaXtrain, sigmaYtrain=sigmaYtrain,
				 X.score=T, X.score.low=T.low, X.loading=P, Y.loading=Q, X.weight=W, 
				 residuals=residuals, residuals.nc=residuals.nc,
				 hatY=hatY, hatY.nc=hatY.nc,
				 hatYtest=hatYtest, hatYtest.nc=hatYtest.nc,
				 A=A, betamat=betamat, new2As=new2As,
				 lambda.l1=lambda.l1, ncomp=ncomp,
				 V=V, adapt=adapt)

\value{
	A list with the following components:
	\item{Xtrain}{the design matrix.}
	\item{Ytrain}{the response observations.}
	\item{sXtrain}{the centered if so and scaled if so design matrix.}
	\item{sYtrain}{the centered if so and scaled if so response.}
	\item{betahat}{the linear coefficients in model \code{sYtrain = sXtrain \%*\% betahat + residuals}.}
	\item{betahat.nc}{the (p+1) vector containing the coefficients and intercept for the non centered
	non scaled model \code{Ytrain = cbnd(rep(1,ntrain),Xtrain) \%*\% betahat.nc + residuals.nc}.}
	\item{meanXtrain}{the (p) vector of Xtrain column mean, used for centering if so.}
	\item{sigmaXtrain}{the (p) vector of Xtrain column standard deviation, used for scaling if so.}
	\item{meanYtrain}{the mean of Ytrain, used for centering if so.}
	\item{sigmaYtrain}{the standard deviation of Ytrain, used for centering if so.}
	\item{X.score}{a (n x ncomp) matrix being the observations coordinates or scores in the 
	new component basis produced by the compression step (sparse PLS). Each column t_k of 
	\code{X.score} is a new component.}
	\item{X.weight}{a (p x ncomp) matrix being the coefficients of predictors in each components 
	produced by sparse PLS. Each column w_k of \code{X.weight} verifies t_k = Xtrain x w_k (as a matrix product).}
	
	\item{hatY}{the (ntrain) vector containing the estimated reponse value on the train 
	set of predictors Xtrain.}
	\item{hatYtest}{the (ntest) vector containing the predicted labels for the observations from 
	\code{Xtest} if non null.}
	\item{DeletedCol}{the vector containing the column number of \code{Xtrain} when the 
	variance of the corresponding predictor variable is null. Otherwise \code{DeletedCol}=NULL}
	\item{A}{the active set of predictors selected by the procedures. \code{A} is a subset of 1:p}
	\item{converged}{a \{0,1\} value indicating whether the IRLS algorithm converged 
	in less than \code{maxIter} iterations or not.}


	\item{lambda.ridge}{the ridge hyper-parameter used to fit the model.}
	\item{lambda.l1}{the sparse hyper-parameter used to fit the model.}
	\item{ncomp}{the number of components used to fit the model.}
	\item{V}{the (ntrain x ntrain) matrix used to weight the metric in the sparse PLS step.\code{V} 
	is the inverse of the covariance matrix of the pseudo-response produced by the RIRLS step.}
	\item{proba.test}{the (ntest) vector of estimated probabilities for the observations in \code{Xtest}, 
	used to predict the \code{hatYtest} labels.}
}

\references{
 
	G. Durif, F. Picard, S. Lambert-Lacroix (2015). Adaptive sparse PLS for logistic regression, (in prep), 
	available on (\url{http://arxiv.org/}). 

}

\author{
	Ghislain Durif (\url{http://lbbe.univ-lyon1.fr/-Durif-Ghislain-.html}). 

}
\seealso{\code{\link{rirls.spls.tune}}.}

\examples{
### load plsgenomics library
library(plsgenomics)

### generating data
n = 100
p = 1000
sample1 = sample.bin(n=100, p=1000, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, mean.H=0.2, sigma.H=10, mean.F=0, sigma.F=5)

X = sample1$X
Y = sample1$Y

### splitting between learning and testing set
index.train = sort(sample(1:n, size=round(0.7*n)))
index.test = (1:n)[-index.train]

Xtrain = X[index.train,]
Ytrain = Y[index.train,]

Xtest = X[index.test,]
Ytest = Y[index.test,]

### fitting the model, and predicting new observations
model1 = rirls.spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.ridge=2, lambda.l1=0.5, ncomp=2, Xtest=Xtest, adapt=TRUE, maxIter=100, svd.decompose=TRUE)
str(model1)

### prediction error rate
sum(model1$hatYtest!=Ytest) / length(index.test)

}
