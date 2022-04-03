\name{poLCA.simdata}
\alias{poLCA.simdata}
\title{Create simulated cross-classification data}
\description{Uses the latent class model's assumed data-generating process to create a simulated dataset that can be used to test the properties of the poLCA latent class and latent class regression estimator.}
\usage{
poLCA.simdata(N = 5000, probs = NULL, nclass = 2, ndv = 4, 
              nresp = NULL, x = NULL, niv = 0, b = NULL, 
              P = NULL, missval = FALSE, pctmiss = NULL)
}
\arguments{
\item{N}{number of observations.}
\item{probs}{a list of matrices of dimension \code{nclass} by \code{nresp} with each matrix corresponding to one manifest variable, and each row containing the class-conditional outcome probabilities (which must sum to 1)  If \code{probs} is \code{NULL} (default) then the outcome probabilities are generated randomly.}
\item{nclass}{number of latent classes. If\code{probs} is specified, then \code{nclass} is set equal to the number of rows in each matrix in that list. If \code{P} is specified, then \code{nclass} is set equal to the length of that vector. If \code{b} is specified, then \code{nclass} is set equal to one greater than the number of columns in \code{b}. Otherwise, the default is two.}
\item{ndv}{number of manifest variables.  If \code{probs} is specified, then \code{ndv} is set equal to the number of matrices in that list. If \code{nresp} is specified, then \code{ndv} is set equal to the length of that vector. Otherwise, the default is four.}
\item{nresp}{number of possible outcomes for each manifest variable. If \code{probs} is specified, then \code{ndv} is set equal to the number of columns in each matrix in that list. If both \code{probs} and \code{nresp} are \code{NULL} (default), then the manifest variables are assigned a random number of outcomes between two and five.}
\item{x}{a matrix of concomicant variables with \code{N} rows and \code{niv} columns. If \code{x=NULL} (default), but \code{niv>0}, then \code{niv} concomitant variables will be  generated as mutually independent random draws from a standard normal distribution.}
\item{niv}{number of concomitant variables (covariates).  Setting \code{niv=0} (default) creates a data set assuming no covariates.  If \code{nclass=1} then \code{niv} is automatically set equal to 0. If both \code{x} and \code{niv} are entered, then the number of columns in \code{x} overrides the value of \code{niv}.  The number of rows in \code{b}, less one, also overrides \code{niv}.}
\item{b}{when using covariates, an \code{niv+1} by \code{nclass-1} matrix of (multinomial) logit coefficients. If \code{b} is \code{NULL} (default), then coefficients are generated as random integers between -2 and 2.}
\item{P}{a vector of mixing proportions (class population shares) of length \code{nclass}. \code{P} must sum to 1.  Disregarded if \code{b} is specified or \code{niv>1} because then \code{P} is, in part, a function of the concomitant variables.  If \code{P} is \code{NULL} (default), then the mixing proportions are generated randomly.}
\item{missval}{logical. If \code{TRUE} then a fraction \code{pctmiss} of the manifest variables are randomly dropped as missing values.  Default is \code{FALSE}.}
\item{pctmiss}{percentage of values to be dropped as missing, if \code{missval=TRUE}.  If \code{pctmiss} is \code{NULL} (default), then a value between 5 and 40 percent is chosen randomly.}
}
\details{Note that entering \code{probs} overrides \code{nclass}, \code{ndv}, and \code{nresp}.  It also overrides \code{P} if the length of the \code{P} vector is not equal to the length of the \code{probs} list.  Likewise, if \code{probs=NULL}, then \code{length(nresp)} overrides \code{ndv} and \code{length(P)} overrides \code{nclass}.  Setting \code{niv>1} causes any user-entered value of \code{P} to be disregarded.}
\value{
\item{dat}{a data frame containing the simulated variables.  Variable names for manifest variables are Y1, Y2, etc.  Variable names for concomitant variables are X1, X2, etc.}
\item{probs}{a list of matrices of dimension \code{nclass} by \code{nresp} containing the class-conditional response probabilities.}
\item{nresp}{a vector containing the number of possible outcomes for each manifest variable.}
\item{b}{coefficients on covariates, if used.}
\item{P}{mixing proportions corresponding to each latent class.}
\item{pctmiss}{percent of observations missing.}
\item{trueclass}{\code{N} by 1 vector containing the "true" class membership for each individual.}
}
\seealso{\code{\link{poLCA}}}
\examples{
# Create a sample data set with 3 classes and no covariates 
# and run poLCA to recover the specified parameters.

probs <- list(matrix(c(0.6, 0.1, 0.3,
                       0.6, 0.3, 0.1,
		       0.3, 0.1, 0.6),
		     ncol = 3,byrow = TRUE), # conditional resp prob to Y1
              matrix(c(0.2, 0.8,
	               0.7, 0.3,
		       0.3, 0.7),
		     ncol = 2, byrow = TRUE), # conditional resp prob to Y2
              matrix(c(0.3, 0.6, 0.1,
	               0.1, 0.3, 0.6,
		       0.3, 0.6, 0.1),
		     ncol = 3,byrow = TRUE), # conditional resp prob to Y3
              matrix(c(0.1, 0.1, 0.5, 0.3,
	               0.5, 0.3, 0.1, 0.1,
		       0.3, 0.1, 0.1, 0.5),
		     ncol = 4,byrow = TRUE), # conditional resp prob to Y4
              matrix(c(0.1, 0.1, 0.8,
	               0.1, 0.8, 0.1,
		       0.8, 0.1, 0.1),
		     ncol = 3,
		     byrow = TRUE)) # conditional resp prob to Y5
simdat <- poLCA.simdata(N=1000,probs,P=c(0.2,0.3,0.5))
f1 <- cbind(Y1,Y2,Y3,Y4,Y5)~1
lc1 <- poLCA(f1,simdat$dat,nclass=3)
table(lc1$predclass,simdat$trueclass)

# Create a sample dataset with 2 classes and three covariates.
# Then compare predicted class memberships when the model is 
# estimated "correctly" with covariates to when it is estimated
# "incorrectly" without covariates.
simdat2 <- poLCA.simdata(N=1000,ndv=7,niv=3,nclass=2,b=matrix(c(1,-2,1,-1)))
f2a <- cbind(Y1,Y2,Y3,Y4,Y5,Y6,Y7)~X1+X2+X3
lc2a <- poLCA(f2a,simdat2$dat,nclass=2)
f2b <- cbind(Y1,Y2,Y3,Y4,Y5,Y6,Y7)~1
lc2b <- poLCA(f2b,simdat2$dat,nclass=2)
table(lc2a$predclass,lc2b$predclass)
}
\keyword{methods}
