\name{carcinoma}
\alias{carcinoma}
\docType{data}
\title{Diagnoses of carcinoma (sample data)}
\description{Dichotomous ratings by seven pathologists of 118 slides for the presence or absence of carcinoma in the uterine cervix. Pathologists are labeled \code{A} through \code{G}. There were 20 different observed response patterns. This data set appears in Agresti (2002, p. 542) as Table 13.1.}
\usage{data(carcinoma)}
\format{A data frame with 118 observations on 7 variables representing pathologist ratings with 1 denoting "no" and 2 denoting "yes".}
\source{Agresti, Alan. 2002. \emph{Categorical Data Analysis, second edition}. Hoboken: John Wiley \& Sons.}
\examples{
##
## Replication of latent class models in Agresti (2002, p. 543), 
## Table 13.2 and Table 13.3.
##
data(carcinoma)
f <- cbind(A,B,C,D,E,F,G)~1
lca2 <- poLCA(f,carcinoma,nclass=2) # log-likelihood: -317.2568
lca3 <- poLCA(f,carcinoma,nclass=3) # log-likelihood: -293.705
lca4 <- poLCA(f,carcinoma,nclass=4,nrep=10,maxiter=5000) # log-likelihood: -289.2858 
}
\keyword{datasets}
