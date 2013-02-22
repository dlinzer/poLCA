\name{poLCA.predcell}
\alias{poLCA.predcell}
\title{Predicted cell percentages in a latent class model}
\description{Calculates the predicted cell percentages from a latent class model, for specified values of the manifest variables.}
\usage{ poLCA.predcell(lc,y) }
\arguments{
  \item{lc}{A model object estimated using the \code{poLCA} function.}
  \item{y}{A vector or matrix containing series of responses on the manifest variables in \code{lc}.}
}
\details{
The parameters estimated by a latent class model can be used to produce a density estimate of the underlying probability mass function across the cells in the multi-way table of manifest variables.  This function calculates cell percentages for that density estimate, corresponding to selected sets of responses on the manifest variables, \code{y}.
}
\value{
A vector containing cell percentages corresponding to the specified sets of responses \code{y}, based on the estimated latent class model \code{lc}.
}
\seealso{ \code{\link{poLCA}} }
\examples{
data(carcinoma)
f <- cbind(A,B,C,D,E,F,G)~1
lca3 <- poLCA(f,carcinoma,nclass=3) # log-likelihood: -293.705

# Only 20 out of 32 possible response patterns are observed
lca3$predcell

# Produce cell probabilities for one sequence of responses
poLCA.predcell(lc=lca3,y=c(1,1,1,1,1,1,1))

# Estimated probabilities for a cell with zero observations
poLCA.predcell(lc=lca3,y=c(1,1,1,1,1,1,2))

# Cell probabilities for both cells at once; y entered as a matrix
poLCA.predcell(lc=lca3,y=rbind(c(1,1,1,1,1,1,1),c(1,1,1,1,1,1,2)))
}
\keyword{methods}
