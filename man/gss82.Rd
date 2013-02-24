\name{gss82}
\alias{gss82}
\docType{data}
\title{1982 General Social Survey (sample data)}
\description{
Attitudes towards survey taking across two dichotomous and two trichotomous items among 1202 white respondents to the 1982 General Social Survey. Respondents give their opinion of the purpose of surveys (\code{PURPOSE}; good/depends/waste of time and money), the accuracy of surveys (\code{ACCURACY}; mostly true/not true), their understanding of survey questions (\code{UNDERSTA}; good/fair, poor), and how well they cooperated with the interviewer (\code{COOPERAT}; interested/cooperative/impatient, hostile). This data set appears in McCutcheon (1987, p. 30) as Table 3.1.
}
\usage{data(gss82)}
\format{A data frame with 1202 observations on 4 survey variables.}
\source{McCutcheon, A.L. 1987. \emph{Latent class analysis}. Newbury Park: SAGE Publications.}
\examples{
data(gss82)
f <- cbind(PURPOSE,ACCURACY,UNDERSTA,COOPERAT)~1
gss.lc2 <- poLCA(f,gss82,nclass=2) # log-likelihood = -2783.268

# Could also try:
# gss.lc3 <- poLCA(f,gss82,nclass=3,maxiter=3000,nrep=10) # log-likelihood = -2754.545
# gss.lc4 <- poLCA(f,gss82,nclass=4,maxiter=15000,nrep=10,tol=1e-7) # log-likelihood = -2746.621
}
\keyword{datasets}
