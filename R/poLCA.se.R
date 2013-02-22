poLCA.se <-
function(y,x,probs,prior,rgivy) {
    J <- ncol(y)
    R <- ncol(prior)
    K.j <- sapply(probs,ncol)
    N <- nrow(y)
    ymat <- y
    y <- list()
    for (j in 1:J) {
        y[[j]] <- matrix(0,nrow=N,ncol=K.j[j])
        y[[j]][cbind(c(1:N),ymat[,j])] <- 1
        y[[j]][ymat[,j]==0,] <- NA      # To handle missing values
    }
    s <- NULL
    # score matrix contains \sum_j [R(K.j-1)] columns correpsonding to log-odds response probs...
    for (r in 1:R) {
        for (j in 1:J) {
            s <- cbind(s,rgivy[,r] * t(t(y[[j]][,2:K.j[j]]) - probs[[j]][r,2:K.j[j]]))
        }
    }
    # ...and (R-1)*ncol(x) columns corresponding to coefficients (betas)
    ppdiff <- rgivy-prior
    if (R>1) for (r in 2:R) { s <- cbind(s,x*ppdiff[,r]) }
    
    s[is.na(s)] <- 0      # To handle missing values
    info <- t(s) %*% s    # Information matrix
    VCE <- ginv(info)     # VCE matrix of log-odds response probs and covariate coefficients

    # Variance of class conditional response probs using delta fn. transformation with 
    # Jacobian a block diagonal matrix (across r) of block diagonal matrices (across j)
    VCE.lo <- VCE[1:sum(R*(K.j-1)),1:sum(R*(K.j-1))]
    Jac <- matrix(0,nrow=nrow(VCE.lo)+(J*R),ncol=ncol(VCE.lo))
    rpos <- cpos <- 1
    for (r in 1:R) {
        for (j in 1:J) {
            Jsub <- -(probs[[j]][r,] %*% t(probs[[j]][r,]))
            diag(Jsub) <- probs[[j]][r,]*(1-probs[[j]][r,])
            Jsub <- Jsub[,-1]
            Jac[rpos:(rpos+K.j[j]-1),cpos:(cpos+K.j[j]-2)] <- Jsub
            rpos <- rpos+K.j[j]
            cpos <- cpos+K.j[j]-1
        }
    }
    VCE.probs <- Jac %*% VCE.lo %*% t(Jac)
    
    maindiag <- diag(VCE.probs)
    maindiag[maindiag<0] <- 0 # error trap
    se.probs.vec <- sqrt(maindiag)
    se.probs <- list()
    for (j in 1:J) { se.probs[[j]] <- matrix(0,0,K.j[j]) }
    pos <- 1
    for (r in 1:R) {
        for (j in 1:J) { 
            se.probs[[j]] <- rbind(se.probs[[j]],se.probs.vec[pos:(pos+K.j[j]-1)])
            pos <- pos+K.j[j]
        }
    }

    # Variance of mixing proportions (priors) and coefficients (betas)
    if (R>1) {
        VCE.beta <- VCE[(1+sum(R*(K.j-1))):dim(VCE)[1],(1+sum(R*(K.j-1))):dim(VCE)[2]]
        se.beta <- matrix(sqrt(diag(VCE.beta)),nrow=ncol(x),ncol=(R-1))
    
        ptp <- array(NA,dim=c(R,R,N))
        for (n in 1:N) {
            ptp[,,n] <- -(prior[n,] %*% t(prior[n,]))
            diag(ptp[,,n]) <- prior[n,] * (1-prior[n,])
        }
        Jac.mix <- NULL
        for (r in 2:R) {
            for (l in 1:ncol(x)) {
                Jac.mix <- cbind(Jac.mix,colMeans(t(ptp[,r,]) * x[,l]))
            }
        }
        VCE.mix <- Jac.mix %*% VCE.beta %*% t(Jac.mix)
        se.mix <- sqrt(diag(VCE.mix))
    } else {
        VCE.beta <- se.beta <- se.mix <- NULL
    }
    return( list(probs=se.probs,P=se.mix,b=se.beta,var.b=VCE.beta) )
}
