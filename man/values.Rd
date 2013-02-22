\name{values}
\alias{values}
\docType{data}
\title{Universalistic vs. particularistic values (sample data)}
\description{Dichotomous survey responses from 216 respondents to four questions (\code{A}, \code{B}, \code{C}, \code{D}) measuring tendencies towards "universalistic" or "particularistic" values.  This data set appears in Goodman (2002, p. 14) as Table 4, and previously appeared in Goodman (1974) and Stouffer and Toby (1951).}
\usage{data(values)}
\format{A data frame with 216 observations on 4 variables representing survey responses to dichotomous questions, with 1 denoting the "particularistic" values response and 2 denoting the "universalistic" values response.}
\source{
Stouffer, S.A. and J. Toby. 1951. "Role conflict and personality." \emph{American Journal of Sociology}. 56: 395:406.

Goodman, Leo A. 1974. "Exploratory Latent-Structure Analysis Using Both Identifiable and Unidentifiable Models." \emph{Biometrika}. 61(2): 215-231.

Goodman, Leo A. 2002. "Latent Class Analysis; The Empirical Study of Latent Types, Latent Variables, and Latent Structures." in Jacques A. Hagenaars and Allan L. McCutcheon, eds. \emph{Applied Latent Class Analysis}. Cambridge: Cambridge University Press.
}
\examples{
##
## Replication of latent class models in Goodman (2002), 
## Tables 5b, 5c, and 6.
##
data(values)
f <- cbind(A,B,C,D)~1
M0 <- poLCA(f,values,nclass=1) # log-likelihood: -543.6498
M1 <- poLCA(f,values,nclass=2) # log-likelihood: -504.4677
M2 <- poLCA(f,values,nclass=3,maxiter=8000) # log-likelihood: -503.3011
}
\keyword{datasets}
