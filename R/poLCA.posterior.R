poLCA.posterior <-
function(lc,y,x=NULL) {
    if (is.vector(y) | any(dim(y)==1)) {
        y <- matrix(y,nrow=1)
    }
    y[is.na(y)] <- 0

    if (!is.null(x)) {
        if (is.vector(x)) x <- matrix(x,nrow=1)
        x <- cbind(1,x)
    }

    if ((ncol(lc$x)>1) & (!is.null(x))) {
        prior <- poLCA.updatePrior(lc$coeff,x,length(lc$P))
    } else {
        prior <- matrix(lc$P,nrow=nrow(y),ncol=length(lc$P),byrow=T)
    }
    ret <- poLCA.postClass.C(prior,poLCA.vectorize(lc$probs),y)
    return(ret)
}
