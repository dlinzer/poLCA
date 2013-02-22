poLCA.makeplot.dich <-
function(probs,P,y,ti) {
    R <- nrow(probs[[1]])
    pi.class <- matrix(NA,nrow=length(probs),ncol=R)
    for (j in 1:length(probs)) {
        pi.class[j,] <- probs[[j]][,2]
    }
    dimnames(pi.class) <- list(names(y),round(P,4))
    ds.plot <- data.frame(Classes=as.vector(col(pi.class)),Manifest.variables=as.vector(row(pi.class)),value=as.vector(pi.class))
    vis <- scatterplot3d(ds.plot,type="h",lwd=5,pch=" ",x.ticklabs=colnames(pi.class),y.ticklabs=colnames(y),z.ticklabs=" ",
            xlab="Classes; population share",ylab="Manifest variables",zlab="Pr(outcome)",color=2,main=ti,y.margin.add=0.2,
            mar=c(6,3,3,3),lab=c(R-1,ncol(y)-1),zlim=c(0,1),box=FALSE,cex.main=1,angle=83)
}
