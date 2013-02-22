poLCA.postClass.C <-
function(prior,vp,y) {
    ret <-  .C("postclass",
                as.double(t(prior)),
                as.double(vp$vecprobs),
                as.integer(t(y)),
                as.integer(length(vp$numChoices)),
                as.integer(dim(y)[1]),
                as.integer(vp$numChoices),
                as.integer(vp$classes),
                posterior = double(dim(y)[1]*vp$classes)
            )
    ret$posterior <- matrix(ret$posterior,ncol=vp$classes,byrow=TRUE)
    return(ret$posterior)
}
