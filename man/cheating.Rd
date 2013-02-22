\name{cheating}
\alias{cheating}
\docType{data}
\title{GPA and chronic cheating (sample data)}
\description{
Dichotomous responses by 319 undergraduates to four questions about cheating behavior, and each student's academic GPA.  

Students responded either (1) no or (2) yes as to whether they had ever lied to avoid taking an exam (\code{LIEEXAM}), lied to avoid handing a term paper in on time (\code{LIEPAPER}), purchased a term paper to hand in as their own or had obtained a copy of an exam prior to taking the exam (\code{FRAUD}), or copied answers during an exam from someone sitting near to them (\code{COPYEXAM}).

The \code{GPA} variable is partitioned into five groups: (1) 2.99 or less; (2) 3.00-3.25; (3) 3.26-3.50; (4) 3.51-3.75; (5) 3.76-4.00.

This data set appears in Dayton (1998, pp. 33 and 85) as Tables 3.4 and 7.1.
}
\usage{data(cheating)}
\format{A data frame with 319 observations on 5 variables. Note: GPA data were not available for four students who reported never cheating.}
\source{Dayton, C. Mitchell. 1998. \emph{Latent Class Scaling Analysis}. Thousand Oaks, CA: SAGE Publications.}
\examples{
##
## Replication of latent class models in Dayton (1998)
##
## Example 1. Two-class LCA. (Table 3.3, p. 32)
##
data(cheating)
f <- cbind(LIEEXAM,LIEPAPER,FRAUD,COPYEXAM)~1
ch2 <- poLCA(f,cheating,nclass=2)	# log-likelihood: -440.0271 

##
## Example 2. Two-class latent class regression using 
## GPA as a covariate to predict class membership as 
## "cheaters" vs. "non-cheaters".
## (Table 7.1, p. 85, and Figure 7.1, p. 86)
##
f2 <- cbind(LIEEXAM,LIEPAPER,FRAUD,COPYEXAM)~GPA
ch2c <- poLCA(f2,cheating,nclass=2)	# log-likelihood: -429.6384 
GPAmat <- cbind(1,c(1:5))
exb <- exp(GPAmat \%*\% ch2c$coeff)
matplot(c(1:5),cbind(1/(1+exb),exb/(1+exb)),type="l",lwd=2,
	main="GPA as a predictor of persistent cheating",
	xlab="GPA category, low to high",
	ylab="Probability of latent class membership")
text(1.7,0.3,"Cheaters")
text(1.7,0.7,"Non-cheaters")

##
## Compare results from Example 1 to Example 2.
## Non-simultaneous estimation of effect of GPA on latent class
## membership biases the estimated effect in Example 1. 
##
cheatcl <- which.min(ch2$P)
predcc <- sapply(c(1:5),function(v) mean(ch2$posterior[cheating$GPA==v,cheatcl],na.rm=TRUE))
## Having run Ex.2, add to plot:
matplot(c(1:5),cbind(1-predcc,predcc),type="l",lwd=2,add=TRUE)
text(4,0.14,"Cheaters\n (non-simul. estimate)")
text(4,0.87,"Non-cheaters\n (non-simul. estimate)")
}
\keyword{datasets}
