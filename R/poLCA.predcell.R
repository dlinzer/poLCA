poLCA.predcell <-
function(lc,y) {
    K.j <- sapply(lc$probs,ncol)
    trap <- FALSE
    if (is.vector(y) | any(dim(y)==1)) {
        if ((length(y)!=length(K.j)) | (any(y[1:length(K.j)]>K.j)) | (any(y<=0)) | (any(y!=round(y)))) {
            cat("Error: invalid vector (y) of manifest variable values. \n")
            trap <- TRUE
        } else {
            y <- matrix(y,nrow=1)
        }
    } else {
        if ((ncol(y)!=length(K.j)) | (any(apply(y,2,max)>K.j)) | (any(y<=0)) | (any(y!=round(y)))) {
            cat("Error: invalid matrix (y) of manifest variable values. \n")
            trap <- TRUE
        }
    }
    if (trap) {
        invisible(NULL)
    } else {
        ret <- poLCA.ylik.C(poLCA.vectorize(lc$probs),y) %*% lc$P
        return(ret)
    }
}
