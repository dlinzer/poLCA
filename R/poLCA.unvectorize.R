poLCA.unvectorize <-
function(vp) {
    probs <- list()
    idx <- c(0,cumsum(vp$numChoices*vp$classes))
    for (i in 1:length(vp$numChoices)){
        probs[[i]] <- matrix(vp$vecprobs[(idx[i]+1):idx[i+1]],nrow=vp$classes,byrow=TRUE)
    }
    return(probs)
}
