#' Bivariate residuals for latent class models
#' 
#' Calculate the "bivariate residuals" (BVRs) between pairs of variables 
#' in a latent class model.
#' 
#' This function compares the model-implied (expected) counts in the crosstables
#' of all pairs of observed dependent variables to the observed counts. For each
#' pair, it calculates a "chi-square" statistic,
#' 
#' \deqn{\text{BVR} = \sum_{j, j'} \frac{(n_{jj'} - e_{jj'})^2}{e_{jj'}}},
#' 
#' where \eqn{n_{jj'}} are the observed counts for categories \eqn{j} and \eqn{j'} 
#' of the variables being crosstabulated, and \eqn{e_{jj'}} are
#' the expected counts under the latent class model. 
#' 
#' Note that the BVR does not follow an asymptotic chi-square distribution and
#' for accurate p-values, parametric bootstrapping is necessary (Oberski et al. 2013).
#' 
#' @param fit A poLCA fit object
#' @param tol Optional: tolerance for small expected counts
#' @param rescale_to_df Optional: whether to divide the pairwise "chi-square" values by 
#' the degrees of freedom of the local crosstable. Default is TRUE.
#' @return The table of bivariate residuals
#' @author Daniel Oberski (daniel.oberski@gmail.com)
#' @seealso \code{\link{poLCA}} for fitting the latent class model.
#' @references 
#' Oberski, DL, Van Kollenburg, GH and Vermunt, JK (2013). 
#'   A Monte Carlo evaluation of three methods to detect local dependence in binary data latent class models. 
#'   Advances in Data Analysis and Classification 7 (3), 267-279.
#' @examples
#' data(values)
#' f <- cbind(A, B, C, D) ~ 1
#' M0 <- poLCA(f,values, nclass=1, verbose = FALSE) 
#' bvr(M0) # 12.4, 5.7, 8.3, 15.6, ... 
bvr <- function(fit, tol = 1e-3, rescale_to_df = TRUE) {
  stopifnot(class(fit) == "poLCA")
  
  ov_names <- names(fit$predcell)[1:(ncol(fit$predcell) - 2)]
  ov_combn <- combn(ov_names, 2)
  
  get_bvr <- function(ov_pair) {
    form_obs <- as.formula(paste0("observed ~ ", ov_pair[1], " + ", ov_pair[2]))
    form_exp <- as.formula(paste0("expected ~ ", ov_pair[1], " + ", ov_pair[2]))
    
    counts_obs <- xtabs(form_obs, data = fit$predcell)
    counts_exp <- xtabs(form_exp, data = fit$predcell)
    counts_exp <- ifelse(counts_exp < tol, tol, counts_exp) # Prevent Inf/NaN
    
    bvr_df <- prod(dim(counts_exp) - 1)
    bvr_value <- sum((counts_obs - counts_exp)^2 / counts_exp)
    
    if(rescale_to_df) bvr_value <- bvr_value / bvr_df
    
    attr(bvr_value, "df") <- bvr_df
    
    bvr_value
  }
  
  bvr_pairs <- apply(ov_combn, 2, get_bvr)

  attr(bvr_pairs, "rescale_to_df") <- rescale_to_df
  attr(bvr_pairs, "class") <- "dist"
  attr(bvr_pairs, "Size") <- length(ov_names)
  attr(bvr_pairs, "Labels") <- ov_names
  attr(bvr_pairs, "Diag") <- FALSE
  attr(bvr_pairs, "Upper") <- FALSE
  
  bvr_pairs
}
