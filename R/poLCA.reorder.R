poLCA.reorder <-
function(probs,o.new) {
    J <- length(probs)
    for (j in 1:J) probs[[j]] <- probs[[j]][o.new,]
    return(probs)
}
