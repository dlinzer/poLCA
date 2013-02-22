poLCA.probHat.C <-
function(rgivy,y,vp) {
    ret <-  .C("probhat",
                as.integer(t(y)),
                as.double(t(rgivy)),
                as.integer(length(vp$numChoices)),
                as.integer(dim(y)[1]),
                as.integer(vp$numChoices),
                as.integer(vp$classes),
                ph = double(sum(vp$numChoices)*vp$classes)
            )
    return(ret$ph)
}
