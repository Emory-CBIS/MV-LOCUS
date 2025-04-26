#' Calculate BIC for Decomposition Results (Multi-View)
#'
#' @param Y A list of original input matrices, each of size N × p.
#' @param res The output list returned by `multi_view_decomposition()`.
#'
#' @return A numeric vector of BIC values, one for each view.
#' @export
calculate_bic <- function(Y, res) {
  stopifnot(length(Y) == length(res$A))
  
  bic_values <- numeric(length(Y))
  
  for (i in seq_along(Y)) {
    N <- nrow(Y[[i]])
    p <- ncol(Y[[i]])
    
    A <- res$A[[i]]
    S <- res$S[[i]]
    
    residuals <- Y[[i]] - A %*% S
    sigma <- sqrt(1 / (N * p) * sum(residuals^2))
    sigma <- max(sigma, 1e-8)  # Avoid zero variance issues
    
    loglike <- 0
    for (j in 1:N) {
      mean_vec <- as.vector(A[j, ] %*% S)
      loglike <- loglike - 2 * sum(dnorm(Y[[i]][j, ], mean = mean_vec, sd = sigma, log = TRUE))
    }
    
    L0 <- log(N) * sum(abs(S) > 1e-1)
    
    bic_values[i] <- loglike + L0
  }
  
  return(bic_values)
}
