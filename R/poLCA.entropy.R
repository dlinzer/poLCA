# Recommended update to poLCA.entropy function

poLCA.entropy <- 
  function(lc, type = c("poLCA", "Mplus", "LGold", "LGoldR2")) {
    
    type <- match.arg(type)
    
    # Base entropy function
    entropy <- function(x) {-sum(x * log(x), na.rm = TRUE)}
    
    ## "Classic" poLCA entropy ##
    K.j <- sapply(lc$probs,ncol)
    fullcell <- expand.grid(lapply(K.j,seq,from=1))
    P.c <- poLCA.predcell(lc,fullcell)
    
    # "Classic" poLCA
    entropy.polca <- entropy(P.c)
    
    ## Mplus and LatentGold entropy statistics ##
    machine.tolerance <- sqrt(.Machine$double.eps)
    n <- lc$N # number of observations
    k <- length(lc$P) #number of classes
    
    prob.df <- data.frame(lc$posterior) # posterior probabilities of existing combinations of indicators
    prob.unpivot <- stack(prob.df) # unpivot probabilities
    prob <- subset(prob.unpivot, subset = values > machine.tolerance, select = values) # filter near-zero values since Limit{p->0} -p*log(p) = 0
    
    # Mplus
    entropy.mplus <- 1 - (entropy(prob) / (n * log(k)))
    
    # LatentGold
    entropy.lgold <- entropy(prob)
    
    # LatentGold entropy-based R-squared
    error <- entropy.lgold / n # Standardized entropy as "error"
    error.total <- sum(sapply(lc$P, entropy))
    entropy.R2.lgold <- (error.total - error) / error.total
    
    result <- switch(
      type,
      poLCA = entropy.polca,
      Mplus = entropy.mplus,
      LGold = entropy.lgold,
      LGoldR2 = entropy.R2.lgold
    )
    
    return(result)
  }
