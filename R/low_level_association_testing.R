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


prepare_output <- function(permutation_runs, null_dist_fit, p_value, contingency_table, side, n_covariates, precomp_str, B, output_amount) {
  # basic output: z_value, log_fold_change, p_value
  output <- data.frame(z_value = permutation_runs$z_star,
                       log_fold_change = permutation_runs$log_fold_change,
                       p_value = p_value)

  # intermediate output: ks_fit, empirical p-value, contingency table, fit information
  if (output_amount >= 2) {
    ks_fit <- compute_ks_test(z_null = permutation_runs$z_null,
                              dp = null_dist_fit$dp,
                              distribution = "SN")
    p_value_emp <- compute_empirical_p_value(z_star = permutation_runs$z_star,
                                             z_null = permutation_runs$z_null,
                                             side = side)
    output[names(ks_fit)] <- ks_fit
    output[names(null_dist_fit$dp)] <- null_dist_fit$dp
    output[names(contingency_table)] <- contingency_table
    output <- cbind(output, data.frame(p_value_emp = p_value_emp,
                                       n_covariates = n_covariates,
                                       n_iterations = null_dist_fit$n_iterations,
                                       convergence = null_dist_fit$convergence,
                                       B = B,
                                       precomp_summary = factor(precomp_str)))
  }

  # maximum output: resampled test statistics
  if (output_amount >= 3) {
    resampled_stats <- stats::setNames(permutation_runs$z_null,
                                       paste0("z_null_", seq(1, length(permutation_runs$z_null))))

    output <- cbind(output, t(data.frame(resampled_stats)))
    row.names(output) <- NULL
  }

  return(output)
}
