#' Run permutation test
#'
#' Runs the permutation test using the NB score-test-based z statistic.
#'
#' @param expressions a vector of gene expressions
#' @param fitted_means the vector of fitted means, obtained from regressing the expressions onto the technical factors
#' @param ground_truth_treatment_idxs indexes of the "treatment" cells on the ground truth data
#' @param synthetic_treatment_idxs indexes of the "treatment" cells under the resampled data
#' @param gene_theta estimated theta of gene expression distribution
#'
#' @return a list containing (i) z_star (the ground truth test statistic), (ii) z_null (the vector of resampled test statistics), and (iii) log_fold_change (the estimated log fold change)
run_permutation_test <- function(expressions, fitted_means, ground_truth_treatment_idxs, synthetic_treatment_idxs, gene_theta, side, full_output) {
  # compute test statistic on real data
  y <- expressions[ground_truth_treatment_idxs]
  mu <- fitted_means[ground_truth_treatment_idxs]
  log_fold_change <- log(mean(y)) - log(mean(mu))
  z_star <- compute_nb_test_stat(y, mu, gene_theta)
  z_null <- apply(X = synthetic_treatment_idxs, MARGIN = 2, FUN = function(col) {
    compute_nb_test_stat(expressions[col], fitted_means[col], gene_theta)
  })
  return(list(z_star = z_star, z_null = z_null, log_fold_change = log_fold_change))
}


#' Compute NB test statistic
#'
#' @param y vector of expressions
#' @param mu vector of means
#' @param gene_theta estimated theta of gene expression distribution
#'
#' @return the test statistic
compute_nb_test_stat <- function(y, mu, gene_theta) {
  r_mu <- gene_theta * mu
  y_mu <- y * mu
  r_plus_mu <- gene_theta + mu
  sum_y <- sum(y)
  top <- (y_mu + r_mu)/r_plus_mu
  bottom <- r_mu/r_plus_mu
  z <- (sum_y - sum(top))/sqrt(sum(bottom))
  return(z)
}
