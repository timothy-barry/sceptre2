perform_association_test_lowmoi_odm_v2 <- function(mm_odm, grna_group_info, response_grna_group_pairs, B, output_amount, side) {
  # obtain response odm
  response_odm <- mm_odm |> ondisc::get_modality("response")

  # get ordered response-grna pairs, as well as the vector of unique response ids
  response_grna_group_pairs <- response_grna_group_pairs |> dplyr::arrange(response_id, grna_group)
  response_ids <- as.character(unique(response_grna_group_pairs$response_id))
  grna_groups <- as.character(unique(response_grna_group_pairs$grna_group))

  # get the global cell covariates
  global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()

  print("Running gene precomputations.")
  # loop over the response ids, running the precomputation and saving the results
  gene_precomp_list <- lapply(X = response_ids, FUN = function(response_id) {
    # load expression data
    expressions <- as.numeric(response_odm[[response_id,]])

    # perform the response precomputation
    run_response_precomputation_low_level(expressions = expressions[grna_group_info[["grna_specific_idxs"]][["non-targeting"]]],
                                          covariate_matrix = global_cell_covariates[grna_group_info[["grna_specific_idxs"]][["non-targeting"]],])
  }) |> setNames(response_ids)

  # Next, loop over the gRNAs
  res <- lapply(X = grna_groups, FUN = function(grna_group) {
    # initial grna group-specific vectors and values
    n_cells_curr_grna_group <- grna_group_info[["n_cells_per_grna"]][[grna_group]]
    subset_vect <- c(grna_group_info[["grna_specific_idxs"]][[grna_group]],
                     grna_group_info[["grna_specific_idxs"]][["non-targeting"]])
    synthetic_treatment_idxs <- get_grna_permutation_idxs(grna_group_info[["n_cells_per_grna"]],
                                                          grna_group, B)

    # loop over the response ids
    response_ids <- response_grna_group_pairs |>
      dplyr::filter(grna_group == !!grna_group) |>
      dplyr::pull(response_id) |>
      as.character()
    lapply(X = response_ids, FUN = function(response_id) {
      # print message
      print(paste0("Working on gene ", response_id, " and gRNA group ", grna_group, "."))

      # load the expressions (again)
      expressions <- as.numeric(response_odm[[response_id,]])

      # compute the fitted values of the regression
      pieces <- get_pieces_from_response_precomp(response_precomp = gene_precomp_list[[response_id]]$precomp,
                                                 global_cell_covariates = global_cell_covariates)
      response_theta <- pieces$response_theta
      fitted_means <- pieces$fitted_means

      # obtain vectors to pass to permutation test
      curr_expressions <- expressions[subset_vect]
      curr_fitted_means <- fitted_means[subset_vect]
      ground_truth_treatment_idxs <- seq(1, n_cells_curr_grna_group)

      # call the low-level association test function
      perm_runs <- run_permutations(expressions = curr_expressions,
                                    fitted_means = curr_fitted_means,
                                    ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                    synthetic_treatment_idxs = synthetic_treatment_idxs,
                                    response_theta = response_theta)
      # compute the p-value
      p_value <- compute_empirical_p_value(z_star = perm_runs$z_star,
                                           z_null = perm_runs$z_null,
                                           side = side)

      # output result
      out <- data.frame(z_value = perm_runs$z_star,
                        log_fold_change = perm_runs$log_fold_change,
                        p_value = p_value,
                        grna_group = factor(grna_group),
                        response_id = factor(response_id))

      # append additional info if output_amount = 2
      if (output_amount == 2) {
        z_null <- perm_runs$z_null
        z_null_df <- as.data.frame(matrix(data = z_null, nrow = 1))
        colnames(z_null_df) <- paste0("z_null_", seq_along(z_null))
        z_null_df$theta <- response_theta
        z_null_df$precomp_str <- factor(gene_precomp_list[[response_id]]$precomp_str)
        z_null_df[,names(out)] <- out
        out <- z_null_df
      }
      out
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()
  return(res)
}
