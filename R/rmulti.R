rmulti <-
function(p) {
    if (is.vector(p)) p <- t(p)
    n <- nrow(p)
    p <- matrix(p[,-ncol(p)],nrow=n)
    thresh <- matrix(apply(p,1,cumsum),nrow=n,ncol=ncol(p),byrow=TRUE)
    vals <- matrix(runif(n),ncol=ncol(thresh),nrow=n)
    return(rowSums(vals>thresh)+1)
}
