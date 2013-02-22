poLCA.posterior <-
function(lc,y,x=NULL) {
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
    if (!is.null(x)) {
        if (is.vector(x)) x <- matrix(x,nrow=1)
        x <- cbind(1,x)
        if (nrow(y) != nrow(x)) {
            cat("Error: number of rows in x does not match number of rows in y. \n") 
            trap <- TRUE
        }
        if (ncol(x) != ncol(lc$x)) {
            cat("Error: incorrect number of columns in x. \n")
            trap <- TRUE
        }
    }
    if (trap) {
        invisible(NULL)
    } else {
        if ((ncol(lc$x)>1) & (!is.null(x))) {
            prior <- poLCA.updatePrior(lc$coeff,x,length(lc$P))
        } else {
            prior <- matrix(lc$P,nrow=nrow(y),ncol=length(lc$P),byrow=T)
        }
        ret <- poLCA.postClass.C(prior,poLCA.vectorize(lc$probs),y)
        return(ret)
    }
}
