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
\source{The National Election Studies (\url{http://www.electionstudies.org}). THE 2000 NATIONAL ELECTION STUDY [dataset]. Ann Arbor, MI: University of Michigan, Center for Political Studies [producer and distributor].}
\examples{
##
## 2000 American National Election Study analysis
##
## Example 1. Latent class models with one (loglinear independence) to three classes.
##
data(election)
f <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
           MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~1
nes1 <- poLCA(f,election,nclass=1)  # log-likelihood: -18647.31
nes2 <- poLCA(f,election,nclass=2)  # log-likelihood: -17344.92
nes3 <- poLCA(f,election,nclass=3)  # log-likelihood: -16714.66

##
## Example 2. Three-class model with a single covariate (party, age).
##
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

f2b <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
             MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~AGE
nes2b <- poLCA(f2b,election,nclass=3,nrep=5)    # log-likelihood: -16625.96 
agemat <- cbind(1,c(18:90))
exb <- exp(agemat \%*\% nes2b$coeff)
matplot(c(18:90),(cbind(1,exb)/(1+rowSums(exb))),ylim=c(0,1),type="l",
    main="Age as a predictor of candidate affinity class",
    xlab="Age",
    ylab="Probability of latent class membership",lwd=2,col=1)
text(30,0.55,"Other")
text(30,0.33,"Bush affinity")
text(30,0.17,"Gore affinity")

##
## Example 3. Three-class model with covariates 
## age, education, and age*education interaction.
## Graph shows predicted class probabilities by age, 
## for high-school grads versus college grads.
##
f3 <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
            MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~AGE*EDUC
nes3cov <- poLCA(f3,election,nclass=3,nrep=5)   # log-likelihood: -16601.04
predmat <- list()
for (i in 1:3) {
    predmat[[i]] <- matrix(NA,nrow=73,ncol=2)
    ivmat.HS <- cbind(1,c(18:90),3,(c(18:90)*3))
    ivmat.COLL <- cbind(1,c(18:90),6,(c(18:90)*6))
    exb.HS <- exp(ivmat.HS \%*\% nes3cov$coeff)
    exb.COLL <- exp(ivmat.COLL \%*\% nes3cov$coeff)
    predmat[[i]][,1] <- (cbind(1,exb.HS)/(1+rowSums(exb.HS)))[,i]
    predmat[[i]][,2] <- (cbind(1,exb.COLL)/(1+rowSums(exb.COLL)))[,i]
}
matplot(c(18:90),predmat[[1]],ylim=c(0.1,0.6),type="l",lty=c(1,2),col="blue",
    main="Age and Education as predictors of candidate affinity class",
    xlab="Age",
    ylab="Probability of latent class membership",lwd=2)
matplot(c(18:90),predmat[[2]],type="l",lty=c(1,2),col="red",lwd=2,add=TRUE)
matplot(c(18:90),predmat[[3]],type="l",lty=c(1,2),col="limegreen",lwd=2,add=TRUE)
text(27,0.55,"HS: Other"); text(23,0.49,"College:\n Other")
text(27,0.32,"HS: Bush affinity"); text(31,0.27,"College:\n Bush affinity")
text(27,0.15,"HS: Gore affinity"); text(22,0.235,"College:\n Gore affinity")
}
\keyword{datasets}
