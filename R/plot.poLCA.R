plot.poLCA <-
function(x, ...) {
    K.j <- sapply(x$probs,ncol)
    R <- length(x$P)
    if (max(K.j)==2) {
        poLCA.makeplot.dich(x$probs,x$P,x$y,NULL)
    } else {
        layout(matrix(seq(1,(R+1)),R+1,1),heights=c(rep(5,R),1))
        for (r in 1:R) {
            poLCA.makeplot.poly(x$probs,r,x$y,K.j,paste("Class ",r,": population share = ",round(x$P[r],3),sep=""))
        }
    }
    par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
}
