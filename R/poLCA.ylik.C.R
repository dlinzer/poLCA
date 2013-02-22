poLCA.ylik.C <-
function(vp,y) {
    ret <-  .C("ylik",
                as.double(vp$vecprobs),
                as.integer(t(y)),
                as.integer(dim(y)[1]),
                as.integer(length(vp$numChoices)),
                as.integer(vp$numChoices),
                as.integer(vp$classes),
                lik = double(dim(y)[1]*vp$classes)
            )
    ret$lik <- matrix(ret$lik,ncol=vp$classes,byrow=TRUE)
    return(ret$lik)
}
