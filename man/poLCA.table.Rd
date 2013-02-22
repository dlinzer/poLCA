\name{poLCA.table}
\alias{poLCA.table}
\title{Frequency tables of predicted cell counts from latent class analysis}
\description{Calculates predicted cell frequencies based on an estimated latent class model.}
\usage{ poLCA.table(formula, condition, lc) }
\arguments{
  \item{formula}{A formula expression of the form \code{variable ~ 1} for a one-way frequency distribution, or \code{row ~ column} for two way-tables.}
  \item{condition}{A list containing the values of the manifest variables to hold fixed when creating the table specified by the \code{formula} argument. Setting this to an empty list, \code{condition=list()}, conditions on none of the other manifest variables, producing the marginal frequencies.}
  \item{lc}{A model object previously estimated using the \code{poLCA} function.}
}
\details{
This function outputs predicted cell counts for user-specified combinations of the manifest variables, based on a latent class model estimated by the \code{\link{poLCA}} function.  The \code{predcell} table outputted automatically by \code{poLCA} also contains predicted cell frequencies, but only for cells containing at least one observation.  In contrast, \code{poLCA.table} will calculate predicted cell counts for all cells, including those with zero observations.
}
\value{
A vector or table containing the specified frequency distribution.
}
\seealso{ \code{\link{poLCA}} }
\examples{
data(gss82)
f <- cbind(PURPOSE,ACCURACY,UNDERSTA,COOPERAT)~1
gss.lc2 <- poLCA(f,gss82,nclass=2)
gss.lc2$predcell

poLCA.table(formula=COOPERAT~1,condition=list(PURPOSE=3,ACCURACY=1,UNDERSTA=2),lc=gss.lc2)

poLCA.table(formula=COOPERAT~UNDERSTA,condition=list(PURPOSE=3,ACCURACY=1),lc=gss.lc2)

poLCA.table(formula=COOPERAT~UNDERSTA,condition=list(),lc=gss.lc2)
}
\keyword{methods}
