run_permutations_v2 <- function(expressions, mu_hats, ground_truth_treatment_idxs, synthetic_treatment_idxs, response_theta) {
  # compute distillation vectors a and b
  a <- expressions - (expressions * mu_hats + response_theta * mu_hats)/(response_theta + mu_hats)
  b <- (response_theta * mu_hats)/(response_theta + mu_hats)

  # compute the z nulls
  z_null <- compute_resampled_statistics(a, b, synthetic_treatment_idxs)

  # compute z star
  z_star <- compute_resampled_statistics(a, b, matrix(ground_truth_treatment_idxs, ncol = 1))

  # return list
  return(list(z_star = z_star, z_null = z_null))
}
# compute LFC (for now, deactivated)
# y <- expressions[ground_truth_treatment_idxs]
# mu <- fitted_means[ground_truth_treatment_idxs]
# log_fold_change <- log(mean(y)) - log(mean(mu))


run_permutations_v3 <- function(expressions, mu_hats, ground_truth_treatment_idxs, synthetic_treatment_idxs, response_theta, Z) {
  # compute the pieces needed for the score test
  working_resid <- expressions/mu_hats - 1
  w <- mu_hats/(1 + mu_hats/response_theta)

  # compute z nulls
  z_null <- run_glm_perm_score_test_with_ingredients(Z = Z,
                                                     working_resid = working_resid,
                                                     w = w,
                                                     index_mat = synthetic_treatment_idxs)

  # compute z_star
  z_star <- run_glm_perm_score_test_with_ingredients(Z = Z,
                                                     working_resid = working_resid,
                                                     w = w,
                                                     index_mat = matrix(ground_truth_treatment_idxs, ncol = 1))

  return(list(z_star = z_star, z_null = z_null))
}
