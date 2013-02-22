print.poLCA <-
function(x, ...) {
    R <- length(x$P)
    S <- ifelse(is.na(x$coeff[1]),1,nrow(x$coeff))
    cat("Conditional item response (column) probabilities,\n by outcome variable, for each class (row) \n \n")
    print(lapply(x$probs,round,4))
    cat("Estimated class population shares \n", round(x$P,4), "\n \n")
    cat("Predicted class memberships (by modal posterior prob.) \n",round(table(x$predclass)/x$N,4), "\n \n")
    cat("========================================================= \n")
    cat("Fit for", R, "latent classes: \n")
    cat("========================================================= \n")
    if (S>1) {
        for (r in 2:R) {
            cat(r,"/ 1 \n")
            disp <- data.frame(coeff=round(x$coeff[,(r-1)],5),
                               se=round(x$coeff.se[,(r-1)],5),
                               tval=round(x$coeff[,(r-1)]/x$coeff.se[,(r-1)],3),
                               pr=round(1-(2*abs(pt(x$coeff[,(r-1)]/x$coeff.se[,(r-1)],x$resid.df)-0.5)),3))
            colnames(disp) <- c("Coefficient"," Std. error"," t value"," Pr(>|t|)")
            print(disp)
            cat("========================================================= \n")
        }
    }
    cat("number of observations:", x$N, "\n")
    if(x$N != x$Nobs) cat("number of fully observed cases:", x$Nobs, "\n")
    cat("number of estimated parameters:", x$npar, "\n")
    cat("residual degrees of freedom:", x$resid.df, "\n")
    cat("maximum log-likelihood:", x$llik, "\n \n")
    cat("AIC(",R,"): ",x$aic,"\n",sep="")
    cat("BIC(",R,"): ",x$bic,"\n",sep="")
    if (S==1) cat("G^2(",R,"): ",x$Gsq," (Likelihood ratio/deviance statistic) \n",sep="")
    cat("X^2(",R,"): ",x$Chisq," (Chi-square goodness of fit) \n \n",sep="")
    if (x$numiter==x$maxiter) cat("ALERT: iterations finished, MAXIMUM LIKELIHOOD NOT FOUND \n \n")
    if (!x$probs.start.ok) cat("ALERT: error in user-specified starting values; new start values generated \n \n")
    if (x$npar>x$N) cat("ALERT: number of parameters estimated (",x$npar,") exceeds number of observations (",x$N,") \n \n")
    if (x$resid.df<0) cat("ALERT: negative degrees of freedom; respecify model \n \n")
    if (x$eflag) cat("ALERT: estimation algorithm automatically restarted with new initial values \n \n")
    flush.console()
    invisible(x)
}
