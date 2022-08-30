perform_association_test_lowmoi_odm <- function(mm_odm, grna_group_info, response_grna_group_pairs, B, full_output, side) {
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

    x <- sapply(X = grna_groups, FUN = function(grna_group) {
      print(paste0("Working on response ", response_id, " and gRNA group ", grna_group, "."))
      # initial grna group-specific vectors and values
      n_cells_curr_grna_group <- grna_group_info[["n_cells_per_grna"]][[grna_group]]
      subset_vect <- c(grna_group_info[["grna_specific_idxs"]][[grna_group]],
                       grna_group_info[["grna_specific_idxs"]][["non-targeting"]])
      n_cells_curr_de <- n_cells_curr_grna_group + grna_group_info[["n_cells_per_grna"]][["non-targeting"]]

      # generate synthetic indexes
      set.seed(4)
      synthetic_treatment_idxs <- replicate(n = B, expr = sample.int(n = n_cells_curr_de, size = n_cells_curr_grna_group))

      # obtain vectors to pass to permutation test
      curr_expressions <- expressions[subset_vect]
      curr_fitted_means <- fitted_means[subset_vect]
      ground_truth_treatment_idxs <- seq(1, n_cells_curr_grna_group)

      # compute sample sizes in both treatment and control groups
      contingency_table <- get_contingency_table(curr_expressions, ground_truth_treatment_idxs)

      # call the low-level association test function
      perm_out <- run_permutation_test(expressions = curr_expressions,
                                       fitted_means = curr_fitted_means,
                                       ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                       synthetic_treatment_idxs = synthetic_treatment_idxs,
                                       response_theta = response_theta,
                                       side = side,
                                       full_output = full_output)

      prepare_output(permutation_runs = perm_out$permutation_runs,
                     null_dist_fit = perm_out$null_dist_fit,
                     p_value = perm_out$p_value,
                     contingency_table = contingency_table,
                     side = side,
                     precomp_backup = response_precomp$backup,
                     n_covariates = ncol(global_cell_covariates),
                     precomp_str = response_precomp$precomp_str,
                     full_output = 2)

    }) |> t() |> data.table::as.data.table() |> dplyr::mutate(grna_group = grna_groups,
                                                              response_id = response_id)
  }) |> data.table::rbindlist()
}
