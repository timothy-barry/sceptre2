##########################
# Loop over random indices
##########################
perform_association_test_lowmoi_odm_v3 <- function(mm_odm, grna_group_info, response_grna_group_pairs, B, output_amount, side, max_b_per_batch, in_memory) {
  # obtain response odm
  response_odm <- mm_odm |> ondisc::get_modality("response")

  # get ordered response-grna pairs, as well as the vector of unique response ids
  response_grna_group_pairs <- response_grna_group_pairs |> dplyr::arrange(response_id, grna_group)
  response_ids <- as.character(unique(response_grna_group_pairs$response_id))
  grna_groups <- as.character(unique(response_grna_group_pairs$grna_group))

  # get the global cell covariates
  global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()

  # if in_memory, load the relevant submatrix into memory
  if (in_memory) {
    # load matrix
    curr_idxs <- unlist(grna_group_info$grna_specific_idxs)
    gene_exp_mat <- response_odm[[,curr_idxs]] |> as.matrix()
    rownames(gene_exp_mat) <- ondisc::get_feature_ids(response_odm)

    # update idxs
    grna_names <- names(grna_group_info$grna_specific_idxs)
    grna_names <- setNames(grna_names, grna_names)
    idx_end <- cumsum(grna_group_info$n_cells_per_grna)
    idx_start <- c(0, idx_end[-length(idx_end)]) + 1
    names(idx_end) <- names(idx_start) <- grna_names
    new_grna_idxs <- lapply(grna_names, function(grna_name)
      seq(idx_start[grna_name], idx_end[grna_name]))
    grna_group_info$grna_specific_idxs <- new_grna_idxs

    # update global cell covariates
    global_cell_covariates <- global_cell_covariates[curr_idxs,]
  }

  print("Running gene precomputations.")
  # loop over the response ids, running the precomputation and saving the results
  gene_precomp_list <- lapply(X = response_ids, FUN = function(response_id) {
    cat(paste0("Working on gene ", response_id, ".\n"))
    # load expression data
    if (in_memory) {
      expressions <- gene_exp_mat[response_id,]
    } else {
      expressions <- as.numeric(response_odm[[response_id,]])
    }

    # perform the response precomputation
    run_response_precomputation_low_level(expressions = expressions[grna_group_info[["grna_specific_idxs"]][["non-targeting"]]],
                                          covariate_matrix = global_cell_covariates[grna_group_info[["grna_specific_idxs"]][["non-targeting"]],])
  }) |> setNames(response_ids)
  # obtain the batch sizes
  batch_bs <- diff(unique(c(seq(from = 0, to = B, by = max_b_per_batch), B)))

  # next, loop over gRNAs
  res <- lapply(X = grna_groups, FUN = function(grna_group) {
    cat(paste0("Working on gRNA ", grna_group, "\n"))
    # initial grna group-specific vectors and values
    n_cells_curr_grna_group <- grna_group_info[["n_cells_per_grna"]][[grna_group]]
    subset_vect <- c(grna_group_info[["grna_specific_idxs"]][[grna_group]],
                     grna_group_info[["grna_specific_idxs"]][["non-targeting"]])

    # get the response IDs for this gRNA
    response_ids <- response_grna_group_pairs |>
      dplyr::filter(grna_group == !!grna_group) |>
      dplyr::pull(response_id) |> as.character()
    response_ids <- response_ids |> setNames(response_ids)

    set.seed(4)
    # within gRNA, loop over batches
    grna_wise_result <- lapply(seq_along(batch_bs), function(i) {
      curr_B <- batch_bs[i]
      cat(paste0("\tWorking on batch ", i, " of ", length(batch_bs), "\n"))
      curr_synthetic_treatment_idxs <- get_grna_permutation_idxs(grna_group_info[["n_cells_per_grna"]], grna_group, curr_B)
      # loop over the response IDs
      sapply(X = response_ids, FUN = function(response_id) {
        cat(paste0("\t\tWorking on gene ", response_id, "\n"))
        # for the given response id, load the expressions
        if (in_memory) {
          expressions <- gene_exp_mat[response_id,]
        } else {
          expressions <- as.numeric(response_odm[[response_id,]])
        }

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
        perm_runs <- run_permutations_v2(expressions = curr_expressions,
                                         fitted_means = curr_fitted_means,
                                         ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                         synthetic_treatment_idxs = curr_synthetic_treatment_idxs,
                                         response_theta = response_theta)

        # return the left sum
        left_sum <- sum(perm_runs$z_null <= perm_runs$z_star)
        left_sum
      }) |> matrix(nrow = 1) |> as.data.frame()
    }) |> data.table::rbindlist()
    colnames(grna_wise_result) <- response_ids
    p_vals <- compute_empirical_p_value_from_batch_result(grna_wise_result = grna_wise_result,
                                                          B = B, side = side)
    data.frame(p_value = p_vals |> stats::setNames(NULL),
               grna_group = factor(grna_group),
               response_id = response_ids |> stats::setNames(NULL) |> factor())
  }) |> data.table::rbindlist()
  return(res)
}
