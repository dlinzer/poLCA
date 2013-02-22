poLCA.updatePrior <-
function(b,x,R) {
    b <- matrix(b,ncol=(R-1))
    exb <- exp(x %*% b)
    p <- cbind(1,exb)/(rowSums(exb)+1)
    return(p)
}
