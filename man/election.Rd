\name{election}
\alias{election}
\docType{data}
\title{2000 National Election Studies survey (sample data)}
\description{
Survey data from the 2000 American National Election Study.  Two sets of six questions with four responses each, asking respondents' opinions of how well various traits (moral, caring, knowledgable, good leader, dishonest, intelligent) describe presidential candidates Al Gore and George W. Bush.  The responses are (1) Extremely well; (2) Quite well; (3) Not too well; (4) Not well at all.  Many respondents have varying numbers of missing values on these variables.

The data set also includes potential covariates \code{VOTE3}, the respondent's 2000 vote choice (when asked); \code{AGE}, the respondent's age; \code{EDUC}, the respondent's level of education; \code{GENDER}, the respondent's gender; and \code{PARTY}, the respondent's Democratic-Republican partisan identification.

\code{VOTE3} is coded as (1) Gore; (2) Bush; (3) Other.

\code{EDUC} is coded as (1) 8 grades or less; (2) 9-11 grades, no further schooling; (3) High school diploma or equivalency; (4) More than 12 years of schooling, no higher degree; (5) Junior or community college level degree; (6) BA level degrees, no advanced degree; (7) Advanced degree.
   
\code{GENDER} is coded as (1) Male; (2) Female.

\code{PARTY} is coded as (1) Strong Democrat; (2) Weak Democrat; (3) Independent-Democrat; (4) Independent-Independent; (5) Independent-Republican; (6) Weak Republican; (7) Strong Republican.
}
\usage{data(election)}
\format{A data frame with 1785 observations on 17 survey variables.  Of these, 1311 individuals provided responses on all twelve candidate evaluations.}
\source{The National Election Studies (\url{https://electionstudies.org/}). THE 2000 NATIONAL ELECTION STUDY [dataset]. Ann Arbor, MI: University of Michigan, Center for Political Studies [producer and distributor].}
\examples{
# Latent class models with one (loglinear independence) to three classes
data(election)
f <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
           MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~1
nes1 <- poLCA(f,election,nclass=1)  # log-likelihood: -18647.31
nes2 <- poLCA(f,election,nclass=2)  # log-likelihood: -17344.92
nes3 <- poLCA(f,election,nclass=3)  # log-likelihood: -16714.66

# Three-class model with a single covariate (party)
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
\keyword{datasets}
