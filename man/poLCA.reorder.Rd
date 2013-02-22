\name{poLCA.reorder}
\alias{poLCA.reorder}
\title{Reorder latent classes in poLCA}
\description{A helper function to simplify the reordering of latent classes returned by \code{poLCA}.}
\usage{
poLCA.reorder(probs, o.new)
}
\arguments{
  \item{probs}{a list of class-conditional response probabilities previously used as start values to estimate a particular latent class model using \code{poLCA}.}
  \item{o.new}{a vector of length equal to the number of latent classes in \code{probs}, giving the desired reordering of the latent classes.}
}
\details{Because the latent classes outputted by \code{poLCA} are unordered categories, the numerical order of the classes is arbitrary, and is determined solely by the initial values of the EM algorithm. If \code{probs.start} is set to \code{NULL} (the default) when calling \code{poLCA}, then the function generates the starting values randomly in each run, typically rearranging the latent class labels.  The \code{poLCA.reorder} function is a convenient way to manually adjust the order of the latent classes, by changing the order of the \code{probs.start}.  Refitting the latent class model using these reordered start values will produce a model having the desired category labels.}
\value{The function returns a list of matrices containing the rearranged (by row) class-conditional response probabilities.}
\seealso{\code{\link{poLCA}}}
\examples{
##
## Using the "cheating" sample data set, make the larger
## non-cheater class the first ("reference") class in a 
## latent class regression model.  The coefficient on GPA
## now maintains a consistent interpretation.
##
data(cheating)
f2 <- cbind(LIEEXAM,LIEPAPER,FRAUD,COPYEXAM)~GPA
lc.ch <- poLCA(f2,cheating,nclass=2,verbose=FALSE)
probs.start.new <- poLCA.reorder(lc.ch$probs.start,order(lc.ch$P,decreasing=TRUE))
lc.ch <- poLCA(f2,cheating,nclass=2,probs.start=probs.start.new)
}
\keyword{methods}
