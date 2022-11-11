#' Run permutation test
#'
#' @param expressions a vector of response expressions
#' @param fitted_means the vector of fitted means, obtained from regressing the expressions onto the technical factors
#' @param ground_truth_treatment_idxs indexes of the "treatment" cells on the ground truth data
#' @param synthetic_treatment_idxs indexes of the "treatment" cells under the resampled data
#' @param side sidedness of the test (one of "left", "right", or "both")
#' @param full_output return the full output (TRUE) or reduced output (FALSE)?
#' @param response_theta estimated theta of response expression distribution
run_permutation_test <- function(expressions, fitted_means, ground_truth_treatment_idxs, synthetic_treatment_idxs, response_theta, side) {
  permutation_runs <- run_permutations(expressions, fitted_means, ground_truth_treatment_idxs, synthetic_treatment_idxs, response_theta)
  null_dist_fit <- fit_skew_normal(y = permutation_runs$z_null)
  p_value <- compute_skew_normal_p_value(dp = null_dist_fit$dp, z_star = permutation_runs$z_star, side = side)
  out <- return(list(permutation_runs = permutation_runs,
                     null_dist_fit = null_dist_fit,
                     p_value = p_value))
  return(out)
}


#' Run permutations
#'
#' Runs permutations using the NB score-test-based z statistic.
#'
#' @inheritParams run_permutation_test
#'
#' @return a list containing (i) z_star (the ground truth test statistic), (ii) z_null (the vector of resampled test statistics), and (iii) log_fold_change (the estimated log fold change)
run_permutations <- function(expressions, fitted_means, ground_truth_treatment_idxs, synthetic_treatment_idxs, response_theta) {
  # compute test statistic on real data
  y <- expressions[ground_truth_treatment_idxs]
  mu <- fitted_means[ground_truth_treatment_idxs]
  log_fold_change <- log(mean(y)) - log(mean(mu))
  z_star <- compute_nb_test_stat(y, mu, response_theta)
  z_null <- apply(X = synthetic_treatment_idxs, MARGIN = 2, FUN = function(col) {
    compute_nb_test_stat(expressions[col], fitted_means[col], response_theta)
  })
  return(list(z_star = z_star, z_null = z_null, log_fold_change = log_fold_change))
}


#' Compute NB test statistic
#'
#' @param y vector of expressions
#' @param mu vector of means
#' @param response_theta estimated theta of response expression distribution
#'
#' @return the test statistic
compute_nb_test_stat <- function(y, mu, response_theta) {
  r_mu <- response_theta * mu
  y_mu <- y * mu
  r_plus_mu <- response_theta + mu
  sum_y <- sum(y)
  top <- (y_mu + r_mu)/r_plus_mu
  bottom <- r_mu/r_plus_mu
  z <- (sum_y - sum(top))/sqrt(sum(bottom))
  return(z)
}


run_permutations_v2 <- function(expressions, fitted_means, ground_truth_treatment_idxs, synthetic_treatment_idxs, response_theta) {
  # compute z star
  a <- expressions - (expressions * fitted_means + response_theta * fitted_means)/(response_theta + fitted_means)
  b <- (response_theta * fitted_means)/(response_theta + fitted_means)
  z_star <- sum(a[ground_truth_treatment_idxs])/sqrt(sum(b[ground_truth_treatment_idxs]))

  # compute LFC
  y <- expressions[ground_truth_treatment_idxs]
  mu <- fitted_means[ground_truth_treatment_idxs]
  log_fold_change <- log(mean(y)) - log(mean(mu))

  # compute the z nulls
  z_null <- compute_resampled_statistics(a, b, synthetic_treatment_idxs)

  # return list
  return(list(z_star = z_star, z_null = z_null, log_fold_change = log_fold_change))
}
