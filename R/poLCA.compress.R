poLCA.compress <-
function(y) {
    ym.sorted <- y[do.call(order,data.frame(y)),]
    vars <- ncol(ym.sorted)
    datamat <- ym.sorted[1,]
    freq <- 1
    curpos <- 1
    for (i in 2:nrow(ym.sorted)) {
        if (sum(ym.sorted[i,] == ym.sorted[i-1,])==vars) {
            freq[curpos] <- freq[curpos]+1
        } else {
            datamat <- rbind(datamat,ym.sorted[i,])
            freq <- c(freq,1)
            curpos <- curpos+1
        }
    }
    rownames(datamat) <- c(1:length(freq))
    ret <- list(datamat=datamat,freq=freq)
    return(ret)
}
