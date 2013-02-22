\name{rmulti}
\alias{rmulti}
\title{Random draws from a multinomial distribution}
\description{One random draw from a multinomial distribution or list of multinomial distributions.}
\usage{rmulti(p)}
\arguments{\item{p}{matrix of dimension \code{n} by \code{r} containing probabilities, for each row, of drawing each of \code{r} outcomes.  \code{p} may also be entered as a vector, in which case \code{rmulti} treats it as a matrix of dimension \code{n=1} by \code{r}.}}
\value{Returns a vector of length \code{n}.  Each item represents one draw from the multinomial distribution parameterized by the outcome probabilities in each row of \code{p}.}
\note{Each row of matrix \code{p} must sum to 1 or \code{rmulti} will not work properly.}
\examples{
##
## One draw from a three-category multinomial distribution.
##
p1 <- c(0.7,0.2,0.1)
rmulti(p1)

##
## 10,000 draws from a three-category multinomial distribution.
##
n <- 10000
p2 <- matrix(p1,nrow=n,ncol=length(p1),byrow=TRUE)
rmdraws <- rmulti(p2)
table(rmdraws)/n    # should be approximately 0.7, 0.2, 0.1

##
## 10,000 draws from a mixture of three groups of a
## four-category multinomial distribution.
##
group.p <- matrix(c(0.5,0.3,0.2),nrow=n,ncol=3,byrow=TRUE)
group <- rmulti(group.p)
p3 <- t(matrix(NA,nrow=n,ncol=4))
p3[,group==1] <- c(0.7,0.1,0.1,0.1)
p3[,group==2] <- c(0.1,0.7,0.1,0.1)
p3[,group==3] <- c(0.1,0.1,0.1,0.7)
p3 <- t(p3)
rmdraws3 <- rmulti(p3)
table(group,rmdraws3)
table(group,rmdraws3)/rowSums(table(group,rmdraws3))
}  
\keyword{methods}
