#' Run GLM permutation score test
#'
#' @param fit a fitted GLM object
#' @param index_mat a k-by-B matrix of indexes; the indices should be 1-based
#'
#' @return a vector of p-values
#' @examples
#' n <- 5000
#' B <- 10000
#' theta <- 15
#' beta_v <- c(1, 2, 1)
#' Z <- matrix(data = c(rep(1, n),
#'                   x1 = rbinom(n, 1, 0.4),
#'                   x2 = rnorm(n)),
#'            ncol = 3)
#' lin_pred <- as.numeric(Z %*% beta_v)
#' mus <- exp(lin_pred)
#' y <- sapply(mus, function(curr_mu) MASS::rnegbin(n = 1, mu = curr_mu, theta = theta))
#' fit <- glm(y ~ Z + 0, family = MASS::neg.bin(theta))
#' index_mat <- replicate(n = B,
#' expr = sample.int(n = length(y), size = 250))
#' z_scores <- run_glm_perm_score_test(fit, index_mat)
#' hist(z_scores, freq = FALSE, ylim  = c(0, dnorm(0)))
#' xgrid <- seq(-4, 4, 0.01)
#' d_gaus <- dnorm(x = xgrid)
#' lines(xgrid, d_gaus, col = "red")
#' ks.test(z_scores, pnorm)
run_glm_perm_score_test <- function(fit, index_mat) {
  Z <- fit$model[,-1]
  working_resid <- as.numeric(fit$residuals)
  w <- as.numeric(fit$weights)
  run_glm_perm_score_test_with_ingredients(Z, working_resid, w, index_mat - 1L)
}


run_glm_perm_score_test_with_ingredients <- function(Z, working_resid, w, index_mat) {
  # compute Z^T w (to be used throughout)
  ZtW <- sapply(X = seq(1, length(w)), FUN = function(i) w[i] * Z[i,])

  # compute the precision matrix P = Z^t W Z
  P <- ZtW %*% Z

  # compute the spectral decomposition of P
  P_decomp <- eigen(P)

  # obtain U and Lambda^(-1/2)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1/sqrt(P_decomp$values)

  # compute the matrix B = Lambda^(1/2) U^t (Z^t W)
  B <- (Lambda_minus_half * t(U)) %*% ZtW

  # next, compute the vector W M (Y - mu_hat)
  a <- w * working_resid

  # compute the vector of z-scores
  z_scores <- low_level_score_test_vectorized(a = a, B = B, w = w, index_mat = index_mat)
  return(z_scores)
}


# alternate way of computing B using cholesky decomposition; there does not seem to be a difference
# R <- chol(P)
# R_inv_t <- t(backsolve(R, diag(ncol(R))))
# B <- R_inv_t %*% ZtW
