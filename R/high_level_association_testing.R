perform_association_test_lowmoi_odm <- function(mm_odm, grna_group_info, response_grna_group_pairs, B, output_amount, side, sn_approx) {
  # obtain response odm
  response_odm <- mm_odm |> ondisc::get_modality("response")

  # get ordered response-grna pairs, as well as the vector of unique response ids
  response_grna_group_pairs <- response_grna_group_pairs |> dplyr::arrange(response_id, grna_group)
  response_ids <- as.character(unique(response_grna_group_pairs$response_id))

  # get the global cell covariates
  global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()

  # loop over the response ids
  out <- lapply(X = response_ids, FUN = function(response_id) {
    # load expression data
    expressions <- as.numeric(response_odm[[response_id,]])

    # perform the response precomputation
    response_precomp <- run_response_precomputation_low_level(expressions = expressions[grna_group_info[["grna_specific_idxs"]][["non-targeting"]]],
                                                              covariate_matrix = global_cell_covariates[grna_group_info[["grna_specific_idxs"]][["non-targeting"]],])
    precomp_str <- response_precomp$precomp_str

    # compute the fitted values of the regression
    pieces <- get_pieces_from_response_precomp(response_precomp = response_precomp$precomp,
                                               global_cell_covariates = global_cell_covariates)
    response_theta <- pieces$response_theta
    fitted_means <- pieces$fitted_means

    # iterate over the gRNAs for a given response, generating the permutations and fitting the model
    grna_groups <- dplyr::filter(response_grna_group_pairs, response_id == !!response_id) |>
      dplyr::pull(grna_group) |>
      as.character()

    lapply(X = grna_groups, FUN = function(grna_group) {
      print(paste0("Working on response ", response_id, " and gRNA group ", grna_group, "."))
      # initial grna group-specific vectors and values
      n_cells_curr_grna_group <- grna_group_info[["n_cells_per_grna"]][[grna_group]]
      subset_vect <- c(grna_group_info[["grna_specific_idxs"]][[grna_group]],
                       grna_group_info[["grna_specific_idxs"]][["non-targeting"]])
      synthetic_treatment_idxs <- get_grna_permutation_idxs(grna_group_info[["n_cells_per_grna"]],
                                                            grna_group,
                                                            B)
      # obtain vectors to pass to permutation test
      curr_expressions <- expressions[subset_vect]
      curr_fitted_means <- fitted_means[subset_vect]
      ground_truth_treatment_idxs <- seq(1, n_cells_curr_grna_group)

      # compute sample sizes in both treatment and control groups
      contingency_table <- get_contingency_table(curr_expressions, ground_truth_treatment_idxs)

      # call the low-level association test function
      perm_runs <- run_permutations(expressions = curr_expressions,
                                    fitted_means = curr_fitted_means,
                                    ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                    synthetic_treatment_idxs = synthetic_treatment_idxs,
                                    response_theta = response_theta)

      out <- if (sn_approx) {
        null_dist_fit <- fit_skew_normal(y = perm_runs$z_null)
        p_value <- compute_skew_normal_p_value(dp = null_dist_fit$dp,
                                               z_star = perm_runs$z_star, side = side)
        prepare_output(permutation_runs = perm_runs,
                       null_dist_fit = null_dist_fit,
                       p_value = p_value,
                       contingency_table = contingency_table,
                       side = side,
                       n_covariates = ncol(global_cell_covariates),
                       precomp_str = response_precomp$precomp_str,
                       B = B,
                       output_amount = output_amount)
      } else {
        p_value <- compute_empirical_p_value(z_star = perm_runs$z_star,
                                             z_null = perm_runs$z_null,
                                             side = side)
        data.frame(z_value = perm_runs$z_star,
                   log_fold_change = perm_runs$log_fold_change,
                   p_value = p_value)
      }
      out |>
        dplyr::mutate(grna_group = grna_group, response_id = response_id) |>
        data.table::as.data.table()
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()
}
