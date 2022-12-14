% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_test.R
\name{run_glm_perm_score_test}
\alias{run_glm_perm_score_test}
\title{Run GLM permutation score test}
\usage{
run_glm_perm_score_test(fit, index_mat)
}
\arguments{
\item{fit}{a fitted GLM object}

\item{index_mat}{a k-by-B matrix of indexes; the indices should be 0-based}
}
\value{
a vector of p-values
}
\description{
Run GLM permutation score test
}
\examples{
n <- 5000
B <- 10000
theta <- 15
beta_v <- c(1, 2, 1)
Z <- matrix(data = c(rep(1, n),
                  x1 = rbinom(n, 1, 0.4),
                  x2 = rnorm(n)),
           ncol = 3)
lin_pred <- as.numeric(Z \%*\% beta_v)
mus <- exp(lin_pred)
y <- sapply(mus, function(curr_mu) MASS::rnegbin(n = 1, mu = curr_mu, theta = theta))
fit <- glm(y ~ Z + 0, family = MASS::neg.bin(theta))
index_mat <- replicate(n = B,
expr = sample.int(n = length(y), size = 250)) - 1L
z_scores <- run_glm_perm_score_test(fit, index_mat)
hist(z_scores, freq = FALSE, ylim  = c(0, dnorm(0)))
xgrid <- seq(-4, 4, 0.01)
d_gaus <- dnorm(x = xgrid)
lines(xgrid, d_gaus, col = "red")
ks.test(z_scores, pnorm)
}
