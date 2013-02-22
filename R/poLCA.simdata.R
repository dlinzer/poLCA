poLCA.simdata <-
function(N=5000,probs=NULL,nclass=2,ndv=4,nresp=NULL,x=NULL,niv=0,b=NULL,P=NULL,missval=FALSE,pctmiss=NULL) {
    if (is.null(probs)) {
        if (is.null(nresp)) { nresp <- ceiling(runif(ndv,min=1,max=5)) }
        if (!is.null(P))  { nclass <- length(P) }
        if (!is.null(b)) { nclass <- ncol(b)+1 }
        ndv <- length(nresp)
        probs <- list()
        for(i in 1:ndv) {
            probs[[i]] <- matrix(runif(nclass*nresp[i]),nrow=nclass,ncol=nresp[i])
            probs[[i]] <- probs[[i]]/rowSums(probs[[i]])
        }
    } else {
        ndv <- length(probs)
        nclass <- nrow(probs[[1]])
        nresp <- sapply(probs,ncol)
    }
    if (nclass==1) { 
        niv <- 0
        b <- NULL
        P <- 1
        group <- matrix(1,nrow=N,ncol=1)
    } else {
        if (!is.null(x)) {
            niv <- ncol(x)
            if (nrow(x) != N) {
                cat("ALERT: number of rows of x does not equal N; new covariates will be generated randomly. \n \n")
                x <- NULL
            }
            if (ncol(x) != (nrow(b)-1)) {
                cat("ALERT: number of columns of x does not conform to number of rows in b; new covariates will be generated randomly. \n \n")
                x <- NULL
            }
        }
        if (!is.null(b)) { niv <- nrow(b)-1 }
        if (niv > 0) {
            if (is.null(x)) { x <- matrix(rnorm(N*niv),nrow=N,ncol=niv) }
            colnames(x) <- paste("X",c(1:niv),sep="")
            if (is.null(b)) { b <- matrix(round(runif(((nclass-1)*(niv+1)),min=-2,max=2)),nrow=(niv+1)) }
            prior <- poLCA.updatePrior(c(b),cbind(1,x),nclass)
        } else {
            if (nrow(probs[[1]]) != length(P)) {
                P <- runif(nclass)
                P <- P/sum(P)
            }
            prior <- matrix(P,byrow=TRUE,nrow=N,ncol=nclass)
        }
        group <- rmulti(prior)
    }
    y <- rmulti(probs[[1]][group,])
    for (j in 2:ndv) { y <- cbind(y,rmulti(probs[[j]][group,])) }
    colnames(y) <- paste("Y",c(1:ndv),sep="")
    if (niv > 0) { P <- colMeans(poLCA.postClass.C(prior,poLCA.vectorize(probs),y)) }
    if (missval) {
        if (is.null(pctmiss)) pctmiss <- runif(1,min=0.05,max=0.4)
        make.na <- cbind(ceiling(runif(round(pctmiss*N*ndv),min=0,max=N)),ceiling(runif(round(pctmiss*N*ndv),min=0,max=ndv)))
        y[make.na] <- NA
    }
    ret <- list()
    if (is.null(x)) {
        ret$dat <- data.frame(y)
    } else {
        ret$dat <- data.frame(y,x)
    }
    ret$trueclass <- group
    ret$probs <- probs
    ret$nresp <- nresp
    ret$b <- b
    ret$P <- P
    ret$pctmiss <- pctmiss
    return(ret)
}
