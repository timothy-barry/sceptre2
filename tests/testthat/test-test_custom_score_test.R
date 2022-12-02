test_that("multiplication works", {
  n <- 5000
  B <- 10000
  beta_v <- c(1, 2, 1)
  Z <- matrix(data = c(rep(1, n),
                     x1 = rbinom(n, 1, 0.4),
                     x2 = rnorm(n)),
              ncol = 3)
  lin_pred <- as.numeric(Z %*% beta_v)
  mus <- exp(lin_pred)
  y <- sapply(mus, function(curr_mu) rpois(n = 1, lambda = curr_mu))
  fit <- glm(y ~ Z + 0, family = poisson)

  # generate the uncorrelated index mat
  index_mat <- replicate(n = B,
                         expr = sample.int(n = length(y), size = 100)) - 1L
  # generate dense, binary X matrix
  X_mat <- apply(X = index_mat, MARGIN = 2, FUN = function(col) {
    out <- integer(length = n)
    out[col + 1L] <- 1L
    out
  })

  # compute scores
  z_scores_us <- run_glm_perm_score_test(fit, index_mat)
  z_scores_us_2 <- run_glm_perm_score_test(fit, index_mat)
  z_scores_statmod <- statmod::glm.scoretest(fit, X_mat)

  # test equality
  testthat::expect_true(all(abs(z_scores_us - z_scores_statmod) < 1e-5))

  # next, generate the correlated data
  b <- stats::binomial()
  mus <- b$linkinv(lin_pred - 2)
  X <- sapply(X = mus, FUN = function(mu) rbinom(n = 1, size = 1, prob = mu))
  I <- as.matrix(which(X == 1) - 1L, ncol = 1)

  # compute scores
  z_score_us <- run_glm_perm_score_test(fit, I)
  z_score_statmod <- statmod::glm.scoretest(fit, X)

  # test equality
  testthat::expect_true(abs(z_score_us - z_score_statmod ) < 1e-5)
})
