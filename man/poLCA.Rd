\name{poLCA}
\alias{poLCA}
\title{Latent class analysis of polytomous outcome variables}
\description{Estimates latent class and latent class regression models for polytomous outcome variables.}
\usage{
poLCA(formula, data, nclass = 2, maxiter = 1000, graphs = FALSE, 
      tol = 1e-10, na.rm = TRUE, probs.start = NULL, nrep = 1, 
      verbose = TRUE, calc.se = TRUE)
}
\arguments{
\item{formula}{A formula expression of the form \code{response ~ predictors}. The details of model specification are given below.}
\item{data}{A data frame containing variables in \code{formula}. Manifest variables must contain \emph{only} integer values, and must be coded with consecutive values from 1 to the maximum number of outcomes for each variable. All missing values should be entered as \code{NA}.}
\item{nclass}{The number of latent classes to assume in the model. Setting \code{nclass=1} results in \code{poLCA} estimating the loglinear independence model. The default is two.}
\item{maxiter}{The maximum number of iterations through which the estimation algorithm will cycle.}
\item{graphs}{Logical, for whether \code{poLCA} should graphically display the parameter estimates at the completion of the estimation algorithm. The default is \code{FALSE}.}
\item{tol}{A tolerance value for judging when convergence has been reached.  When the one-iteration change in the estimated log-likelihood is less than \code{tol}, the estimation algorithm stops updating and considers the maximum log-likelihood to have been found.}
\item{na.rm}{Logical, for how \code{poLCA} handles cases with missing values on the manifest variables.  If \code{TRUE}, those cases are removed (listwise deleted) before estimating the model. If \code{FALSE}, cases with missing values are retained. Cases with missing covariates are always removed. The default is \code{TRUE}.}
\item{probs.start}{A list of matrices of class-conditional response probabilities to be used as the starting values for the estimation algorithm.  Each matrix in the list corresponds to one manifest variable, with one row for each latent class, and one column for each outcome.  The default is \code{NULL}, producing random starting values.  Note that if \code{nrep>1}, then any user-specified \code{probs.start} values are only used in the first of the \code{nrep} attempts.}
\item{nrep}{Number of times to estimate the model, using different values of \code{probs.start}.  The default is one.  Setting \code{nrep}>1 automates the search for the global---rather than just a local---maximum of the log-likelihood function. \code{poLCA} returns the parameter estimates corresponding to the model with the greatest log-likelihood.}
\item{verbose}{Logical, indicating whether \code{poLCA} should output to the screen the results of the model.  If \code{FALSE}, no output is produced. The default is \code{TRUE}.}
\item{calc.se}{Logical, indicating whether \code{poLCA} should calculate the standard errors of the estimated class-conditional response probabilities and mixing proportions.  The default is \code{TRUE}; can only be set to \code{FALSE} if estimating a basic model with no concomitant variables specified in \code{formula}.}
}
\details{
Latent class analysis, also known as latent structure analysis, is a technique for the analysis of clustering among observations in multi-way tables of qualitative/categorical variables.  The central idea is to fit a model in which any confounding between the manifest variables can be explained by a single unobserved "latent" categorical variable.  \code{poLCA} uses the assumption of local independence to estimate a mixture model of latent multi-way tables, the number of which (\code{nclass}) is specified by the user.  Estimated parameters include the class-conditional response probabilities for each manifest variable, the "mixing" proportions denoting population share of observations corresponding to each latent multi-way table, and coefficients on any class-predictor covariates, if specified in the model.
  
Model specification: Latent class models have more than one manifest variable, so the response variables are \code{cbind(dv1,dv2,dv3...)} where \code{dv#} refer to variable names in the data frame.  For models with no covariates, the formula is \code{cbind(dv1,dv2,dv3)~1}.  For models with covariates, replace the \code{~1} with the desired function of predictors \code{iv1,iv2,iv3...} as, for example, \code{cbind(dv1,dv2,dv3)~iv1+iv2*iv3}.
  
\code{poLCA} treats all manifest variables as qualitative/categorical/nominal -- NOT as ordinal.
}
\value{
\code{poLCA} returns an object of class poLCA; a list containing the following elements:
\item{y}{data frame of manifest variables.}
\item{x}{data frame of covariates, if specified.}
\item{N}{number of cases used in model.}
\item{Nobs}{number of fully observed cases (less than or equal to \code{N}).}
\item{probs}{estimated class-conditional response probabilities.}
\item{probs.se}{standard errors of estimated class-conditional response probabilities, in the same format as \code{probs}.}
\item{P}{sizes of each latent class; equal to the mixing proportions in the basic latent class model, or the mean of the priors in the latent class regression model.}
\item{P.se}{the standard errors of the estimated \code{P}.}
\item{posterior}{matrix of posterior class membership probabilities; also see function \code{\link{poLCA.posterior}}.}
\item{predclass}{vector of predicted class memberships, by modal assignment.}
\item{predcell}{table of observed versus predicted cell counts for cases with no missing values; also see functions \code{\link{poLCA.table}} and \code{\link{poLCA.predcell}}.}
\item{llik}{maximum value of the log-likelihood.}
\item{numiter}{number of iterations until reaching convergence.}
\item{maxiter}{maximum number of iterations through which the estimation algorithm was set to run.}
\item{coeff}{multinomial logit coefficient estimates on covariates (when estimated). \code{coeff} is a matrix with \code{nclass-1} columns, and one row for each covariate.  All logit coefficients are calculated for classes with respect to class 1.}
\item{coeff.se}{standard errors of coefficient estimates on covariates (when estimated), in the same format as \code{coeff}.}
\item{coeff.V}{covariance matrix of coefficient estimates on covariates (when estimated).}
\item{aic}{Akaike Information Criterion.}
\item{bic}{Bayesian Information Criterion.}
\item{Gsq}{Likelihood ratio/deviance statistic.}
\item{Chisq}{Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables.}
\item{time}{length of time it took to run the model.}
\item{npar}{number of degrees of freedom used by the model (estimated parameters).}
\item{resid.df}{number of residual degrees of freedom.}
\item{attempts}{a vector containing the maximum log-likelihood values found in each of the \code{nrep} attempts to fit the model.}
\item{eflag}{Logical, error flag. \code{TRUE} if estimation algorithm needed to automatically restart with new initial parameters. A restart is caused in the event of computational/rounding errors that result in nonsensical parameter estimates.}
\item{probs.start}{A list of matrices containing the class-conditional response probabilities used as starting values in the estimation algorithm. If the algorithm needed to restart (see \code{eflag}), then this contains the starting values used for the final, successful, run.}
\item{probs.start.ok}{Logical. \code{FALSE} if \code{probs.start} was incorrectly specified by the user, otherwise \code{TRUE}.}
\item{call}{function call to \code{poLCA}.}
}
\references{
Agresti, Alan. 2002. \emph{Categorical Data Analysis, second edition}. Hoboken: John Wiley \& Sons.

Bandeen-Roche, Karen, Diana L. Miglioretti, Scott L. Zeger, and Paul J. Rathouz. 1997. "Latent Variable Regression for Multiple Discrete Outcomes." \emph{Journal of the American Statistical Association}. 92(440): 1375-1386.
 
Hagenaars, Jacques A. and Allan L. McCutcheon, eds. 2002. \emph{Applied Latent Class Analysis}. Cambridge: Cambridge University Press.
 
McLachlan, Geoffrey J. and Thriyambakam Krishnan. 1997. \emph{The EM Algorithm and Extensions}. New York: John Wiley \& Sons.
}
\note{
\code{poLCA} uses EM and Newton-Raphson algorithms to maximize the latent class model log-likelihood function. Depending on the starting parameters, this algorithm may only locate a local, rather than global, maximum. This becomes more and more of a problem as \code{nclass} increases. It is therefore highly advisable to run \code{poLCA} multiple times until you are relatively certain that you have located the global maximum log-likelihood. As long as \code{probs.start=NULL}, each function call will use different (random) initial starting parameters.  Alternatively, setting \code{nrep} to a value greater than one enables the user to estimate the latent class model multiple times with a single call to \code{poLCA}, thus conducting the search for the global maximizer automatically.

The term "Latent class regression" (LCR) can have two meanings.  In this package, LCR models refer to latent class models in which the probability of class membership is predicted by one or more covariates.  However, in other contexts, LCR is also used to refer to regression models in which the manifest variable is partitioned into some specified number of latent classes as part of estimating the regression model. It is a way to simultaneously fit more than one regression to the data when the latent data partition is unknown. The \code{\link[flexmix:flexmix]{flexmix}} function in package \pkg{flexmix} will estimate this other type of LCR model.  Because of these terminology issues, the LCR models this package estimates are sometimes termed "latent class models with covariates" or "concomitant-variable latent class analysis," both of which are accurate descriptions of this model.
 
A more detailed presentation is available online at \url{https://www.jstatsoft.org/article/view/v042i10}.
}
\examples{
##
## Three models without covariates:
## M0: Loglinear independence model.
## M1: Two-class latent class model.
## M2: Three-class latent class model.
##
data(values)
f <- cbind(A,B,C,D)~1
M0 <- poLCA(f,values,nclass=1) # log-likelihood: -543.6498
M1 <- poLCA(f,values,nclass=2) # log-likelihood: -504.4677
M2 <- poLCA(f,values,nclass=3,maxiter=8000) # log-likelihood: -503.3011

##
## Three-class model with a single covariate.
##
data(election)
f2a <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
             MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY
nes2a <- poLCA(f2a,election,nclass=3,nrep=5)    # log-likelihood: -16222.32 
pidmat <- cbind(1,c(1:7))
exb <- exp(pidmat \%*\% nes2a$coeff)
matplot(c(1:7),(cbind(1,exb)/(1+rowSums(exb))),ylim=c(0,1),type="l",
    main="Party ID as a predictor of candidate affinity class",
    xlab="Party ID: strong Democratic (1) to strong Republican (7)",
    ylab="Probability of latent class membership",lwd=2,col=1)
text(5.9,0.35,"Other")
text(5.4,0.7,"Bush affinity")
text(1.8,0.6,"Gore affinity")
}
\keyword{methods}
