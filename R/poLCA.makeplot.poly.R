poLCA.makeplot.poly <-
function(probs,r,y,K.j,ti) {
    pi.class <- matrix(NA,nrow=length(probs),ncol=max(K.j))
    for (j in 1:length(probs)) {
        pi.class[j,1:K.j[j]] <- probs[[j]][r,]
    }
    dimnames(pi.class) <- list(as.character(c(1:ncol(y))),as.character(c(1:max(K.j))))
    ds.plot <- data.frame(Manifest.variables=as.vector(row(pi.class)),Outcomes=as.vector(col(pi.class)),value=as.vector(pi.class))
    vis <- scatterplot3d(ds.plot,type="h",lwd=5,pch=" ",x.ticklabs=colnames(y),y.ticklabs=colnames(pi.class),z.ticklabs=" ",
            xlab="Manifest variables",zlab="Pr(outcome)",main=ti,cex.main=1.5,color=2,lab=c(ncol(y)-1,max(K.j)-1),zlim=c(0,1),box=FALSE,
            angle=75,mar=c(3,3,2,3))
}
