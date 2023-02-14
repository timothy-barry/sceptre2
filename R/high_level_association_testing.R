##########################
# Loop over random indices
##########################
perform_association_test_lowmoi_odm_v3 <- function(mm_odm, grna_group_info, response_grna_group_pairs, B, side, max_b_per_batch, in_memory, statistic, return_dist, screen_b) {
  side <- "both"
  lower_thresh <- screen_b * 0.005
  upper_thresh <- 0.995 * screen_b + 1

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
    gene_exp_mat <- response_odm[[response_ids, curr_idxs]] |> as("RsparseMatrix")
    rownames(gene_exp_mat) <- response_ids

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

  cat("Running gene precomputations.\n")
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

  # convert global cell covariates to a matrix
  global_cell_covariates <- as.matrix(global_cell_covariates)
  # obtain the batch sizes
  batch_bs <- diff(unique(c(seq(from = 0, to = B, by = max_b_per_batch), B)))

  # next, loop over gRNAs
  res <- lapply(X = grna_groups, FUN = function(grna_group) {
    cat(paste0("Working on gRNA ", grna_group, "\n"))
    # initial grna group-specific vectors and values
    n_cells_curr_grna_group <- grna_group_info[["n_cells_per_grna"]][[grna_group]]
    subset_vect <- c(grna_group_info[["grna_specific_idxs"]][[grna_group]],
                     grna_group_info[["grna_specific_idxs"]][["non-targeting"]])
    curr_global_cell_covariates <- global_cell_covariates[subset_vect,]
    ground_truth_treatment_idxs <- seq(0, n_cells_curr_grna_group - 1L)

    # get the response IDs for this gRNA
    response_ids <- response_grna_group_pairs |>
      dplyr::filter(grna_group == !!grna_group) |>
      dplyr::pull(response_id) |> as.character()
    response_ids <- response_ids |> setNames(response_ids)
    set.seed(4)

    p_aggregate <- numeric()
    # run the initial screen
    run_screen <- B >= screen_b && !return_dist
    if (run_screen) { # only if B exceeds screen b AND we do not want to return the null distributions
      cat("Running the initial screen.\n")
      screen_res <- perform_association_tests_for_grna(batch_size = screen_b,
                                                       grna_group_info = grna_group_info,
                                                       grna_group = grna_group,
                                                       response_ids = response_ids,
                                                       in_memory = in_memory,
                                                       gene_exp_mat = gene_exp_mat,
                                                       response_odm = response_odm,
                                                       subset_vect = subset_vect,
                                                       curr_global_cell_covariates = curr_global_cell_covariates,
                                                       gene_precomp_list = gene_precomp_list,
                                                       ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                                       statistic = statistic,
                                                       return_dist = return_dist)
      promising_genes_v <- screen_res <= lower_thresh | screen_res >= upper_thresh

      # compute p-values for the filtered genes (if there are any)
      filtered_gene_counts <- screen_res[!promising_genes_v]
      if (length(filtered_gene_counts) >= 1) {
        p_screened <- compute_empirical_p_value_from_batch_result(grna_wise_result = filtered_gene_counts,
                                                                  B = screen_b, side = side)
        p_aggregate <- c(p_aggregate, p_screened)
      }

      # update response_ids to include promising genes only
      promising_genes <- names(which(promising_genes_v))
      response_ids <- promising_genes
    }

    # loop over batches for promising genes
    if (length(response_ids) >= 1) {
      if (run_screen) {
        cat("\tPerforming additional permutations for promising genes.\n")
      } else {
        cat("\tTesting pairs.\n")
      }
      grna_wise_result_promising <- lapply(seq_along(batch_bs), function(i) {
        cat(paste0("\tWorking on batch ", i, " of ", length(batch_bs), "\n"))
        curr_B <- batch_bs[i]
        out <- perform_association_tests_for_grna(batch_size = curr_B,
                                                  grna_group_info = grna_group_info,
                                                  grna_group = grna_group,
                                                  response_ids = response_ids,
                                                  in_memory = in_memory,
                                                  gene_exp_mat = gene_exp_mat,
                                                  response_odm = response_odm,
                                                  subset_vect = subset_vect,
                                                  curr_global_cell_covariates = curr_global_cell_covariates,
                                                  gene_precomp_list = gene_precomp_list,
                                                  ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                                  statistic = statistic,
                                                  return_dist = return_dist) |>
          matrix(ncol = length(response_ids)) |> as.data.frame()
      }) |> data.table::rbindlist()
      if (!return_dist) {
        colnames(grna_wise_result_promising) <- response_ids
        p_promising <- compute_empirical_p_value_from_batch_result(grna_wise_result = grna_wise_result_promising, B = B, side = side)
        p_aggregate <- c(p_aggregate, p_promising)
      }
    }

    # output results in df
    if (!return_dist) {
      out <- data.frame(p_value = p_aggregate |> stats::setNames(NULL),
                        grna_group = factor(grna_group),
                        response_id = names(p_aggregate) |> factor())
    } else {
      grna_wise_result_promising <- data.table::as.data.table(t(grna_wise_result_promising))
      colnames(grna_wise_result_promising) <- c("z_star", paste0("z_", seq(1, B)))
      grna_wise_result_promising$response_id <- factor(response_ids)
      grna_wise_result_promising$grna_group <- factor(grna_group)
      out <- grna_wise_result_promising
    }
    out
  }) |> data.table::rbindlist()
  return(res)
}


perform_association_tests_for_grna <- function(batch_size,
                                               grna_group_info,
                                               grna_group,
                                               response_ids,
                                               in_memory,
                                               gene_exp_mat,
                                               response_odm,
                                               subset_vect,
                                               curr_global_cell_covariates,
                                               gene_precomp_list,
                                               ground_truth_treatment_idxs,
                                               statistic,
                                               return_dist) {
  curr_synthetic_treatment_idxs <- get_grna_permutation_idxs(grna_group_info[["n_cells_per_grna"]], grna_group, batch_size) - 1L
  # loop over the response IDs
  out <- sapply(X = response_ids, FUN = function(response_id) {
    cat(paste0("\t\tWorking on gene ", response_id, "\n"))

    # for the given response id, load the expressions
    expressions <- if (in_memory) gene_exp_mat[response_id, subset_vect] else  as.numeric(response_odm[[response_id, subset_vect]])

    # compute the fitted values of the regression
    pieces <- get_pieces_from_response_precomp(response_precomp = gene_precomp_list[[response_id]]$precomp,
                                               global_cell_covariates = curr_global_cell_covariates)
    response_theta <- pieces$response_theta
    mu_hats <- pieces$mu_hats

    # call the low-level association test function
    if (statistic == "distilled") {
      perm_runs <- run_permutations_v2(expressions = expressions,
                                       mu_hats = mu_hats,
                                       ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                       synthetic_treatment_idxs = curr_synthetic_treatment_idxs,
                                       response_theta = response_theta)
    } else if (statistic == "full") {
      perm_runs <- run_permutations_v3(expressions = expressions,
                                       mu_hats = mu_hats,
                                       ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                       synthetic_treatment_idxs = curr_synthetic_treatment_idxs,
                                       response_theta = response_theta,
                                       Z = curr_global_cell_covariates)
    } else {
      stop("Statistic not recognized.")
    }
    # return the left sum
    if (!return_dist) {
      out <- sum(perm_runs$z_null <= perm_runs$z_star)
    } else {
      out <- c(perm_runs$z_star, perm_runs$z_null)
    }
    return(out)
  })
}
