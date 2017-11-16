poLCA.posterior <-
function(lc,y,x=NULL) {
    if ((ncol(lc$x)>1) & (!is.null(x))) {
        prior <- poLCA.updatePrior(lc$coeff,x,length(lc$P))
    } else {
        prior <- matrix(lc$P,nrow=nrow(y),ncol=length(lc$P),byrow=T)
    }
    ret <- poLCA.postClass.C(prior,poLCA.vectorize(lc$probs),y)
    return(ret)
}
