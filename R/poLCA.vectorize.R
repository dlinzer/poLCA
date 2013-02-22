poLCA.vectorize <-
function(probs) {
    classes <- nrow(probs[[1]])
    vecprobs <- unlist(lapply(probs,t))
    numChoices <- sapply(probs,ncol)
    return(list(vecprobs=vecprobs,numChoices=numChoices,classes=classes))
}
