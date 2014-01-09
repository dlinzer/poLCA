##
## Drew A. Linzer
## drew@votamatic.org
##
## Jeffrey B. Lewis
## jblewis@ucla.edu
##
## January 9, 2014
##
## poLCA: Polytomous variable Latent Class Analysis
##
## Commands to produce examples in user's manual, located at
##
##      http://dlinzer.github.io/poLCA
##

library(poLCA)

#########################################################################
## Section 5.4. Predicted cell frequencies from the latent class model ##
#########################################################################

data(gss82)
f <- cbind(PURPOSE,ACCURACY,UNDERSTA,COOPERAT)~1
gss.lc2 <- poLCA(f,gss82,nclass=2)
gss.lc2$predcell
poLCA.table(formula=COOPERAT~1,condition=list(PURPOSE=3,ACCURACY=1,UNDERSTA=2),lc=gss.lc2)
poLCA.table(formula=COOPERAT~UNDERSTA,condition=list(PURPOSE=3,ACCURACY=1),lc=gss.lc2)


########################################################
## Section 5.7. Recognizing and avoiding local maxima ##
########################################################

data(gss82)
f <- cbind(PURPOSE,ACCURACY,UNDERSTA,COOPERAT)~1

mlmat <- matrix(NA,nrow=500,ncol=4)
for (i in 1:500) { # note, this simulation takes some time to run
    gss.lc <- poLCA(f,gss82,nclass=3,maxiter=3000,tol=1e-7,verbose=F)
    mlmat[i,1] <- gss.lc$llik
    o <- order(gss.lc$probs$UNDERSTA[,1],decreasing=T) # ideal, skeptic, believer
    mlmat[i,-1] <- gss.lc$P[o]
}

# Table 1
tab <- table(round(mlmat[,1],3))

# -2754.545 & & 258 & & 0.621 & 0.172 & 0.207 \\
# -2755.617 & &  14 & & 0.782 & 0.150 & 0.067 \\
# -2755.739 & &  57 & & 0.796 & 0.162 & 0.043 \\
# -2762.005 & &  70 & & 0.508 & 0.392 & 0.099 \\
# -2762.231 & & 101 & & 0.297 & 0.533 & 0.170 \\


######################################################################
## Section 6.1. Basic latent class modeling with the carcinoma data ##
######################################################################

data(carcinoma)
f <- cbind(A,B,C,D,E,F,G)~1
lc2 <- poLCA(f,carcinoma,nclass=2)
lc3 <- poLCA(f,carcinoma,nclass=3)
lc4 <- poLCA(f,carcinoma,nclass=4,maxiter=5000)

print(lc3)

# Figure 1
lc3 <- poLCA(f,carcinoma,nclass=3,graphs=T,maxiter=400)


##########################################################################
## Section 6.2. Latent class regression modeling with the election data ##
##########################################################################

data(election)

## one covariate

f.party <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
                 MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY
nes.party <- poLCA(f.party,election,nclass=3,verbose=F)
# log-likelihood: -16222.32

probs.start <- poLCA.reorder(nes.party$probs.start,order(nes.party$P,decreasing=T))
nes.party <- poLCA(f.party,election,nclass=3,probs.start=probs.start)

# Figure 2
pidmat <- cbind(1,c(1:7)) # matrix of hypothetical party ID values
exb <- exp(pidmat %*% nes.party$coeff)
matplot(c(1:7),(cbind(1,exb)/(1+rowSums(exb))),
                main="Party ID as a predictor of candidate affinity class",
                xlab="Party ID: strong Democratic (1) to strong Republican (7)",
                ylab="Probability of latent class membership",
                ylim=c(0,1),type="l",lwd=3,col=1)
text(5.9,0.35,"Other")
text(5.4,0.7,"Bush affinity")
text(1.8,0.6,"Gore affinity")

## multiple covariates

f.3cov <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
                MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY*AGE
nes.3cov <- poLCA(f.3cov,election,nclass=3,verbose=F)
# log-likelihood: -16135.39 

probs.start <- poLCA.reorder(nes.3cov$probs.start,order(nes.3cov$P,decreasing=T))
nes.3cov <- poLCA(f.3cov,election,nclass=3,probs.start=probs.start)

# Figure 3, left
strdems <- cbind(1,1,c(18:80),(c(18:80)*1))
exb.strdems <- exp(strdems %*% nes.3cov$coeff)
matplot(c(18:80),(cbind(1,exb.strdems)/(1+rowSums(exb.strdems))),
                  main="Age and candidate affinity for strong Democrats",
                  xlab="Age",ylab="Probability of latent class membership",
                  ylim=c(0,1),type="l",col=1,lwd=3)
text(50,0.3,"Other")
text(50,0.05,"Bush affinity")
text(50,0.7,"Gore affinity")

# Figure 3, right
strreps <- cbind(1,7,c(18:80),(c(18:80)*7))
exb.strreps <- exp(strreps %*% nes.3cov$coeff)
matplot(c(18:80),(cbind(1,exb.strreps)/(1+rowSums(exb.strreps))),
                  main="Age and candidate affinity for strong Republicans",
                  xlab="Age",ylab="Probability of latent class membership",
                  ylim=c(0,1),type="l",col=1,lwd=3)
text(50,0.18,"Other")
text(50,0.9,"Bush affinity")
text(50,0.05,"Gore affinity")


# end of file.
