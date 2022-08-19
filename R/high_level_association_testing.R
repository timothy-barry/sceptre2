perform_association_test_lowmoi_odm <- function(mm_odm, grna_group_info, gene_grna_group_pairs, B, full_output, side) {
  # obtain gene odm
  gene_odm <- mm_odm |> ondisc::get_modality("gene")

  # get ordered gene-grna pairs, as well as the vector of unique gene ids
  gene_grna_group_pairs <- gene_grna_group_pairs |> dplyr::arrange(gene_id, grna_group)
  gene_ids <- as.character(unique(gene_grna_group_pairs$gene_id))

  # get the global cell covariates
  global_cell_covariates <- mm_odm |> get_cell_covariates()

  # loop over the gene ids
  out <- lapply(X = gene_ids, FUN = function(gene_id) {
    # load expression data
    expressions <- as.numeric(gene_odm[[gene_id,]])

    # perform the gene precomputation
    gene_precomp <- run_gene_precomputation_low_level(expressions = expressions[grna_group_info[["grna_specific_idxs"]][["non-targeting"]]],
                                                      covariate_matrix = global_cell_covariates[grna_group_info[["grna_specific_idxs"]][["non-targeting"]],])

    # compute the fitted values of the regression
    pieces <- get_pieces_from_gene_precomp(gene_precomp = gene_precomp,
                                           global_cell_covariates = global_cell_covariates)
    gene_theta <- pieces$gene_theta
    fitted_means <- pieces$fitted_means

    # iterate over the gRNAs for a given gene, generating the permutations and fitting the model
    grna_groups <- dplyr::filter(gene_grna_group_pairs, gene_id == !!gene_id) |>
      dplyr::pull(grna_group) |>
      as.character()

    x <- sapply(X = grna_groups, FUN = function(grna_group) {
      print(paste0("Working on gene ", gene_id, " and gRNA group ", grna_group, "."))
      # initial grna group-specific vectors and values
      n_cells_curr_grna_group <- grna_group_info[["n_cells_per_grna"]][[grna_group]]
      subset_vect <- c(grna_group_info[["grna_specific_idxs"]][[grna_group]],
                       grna_group_info[["grna_specific_idxs"]][["non-targeting"]])
      n_cells_curr_de <- n_cells_curr_grna_group + grna_group_info[["n_cells_per_grna"]][["non-targeting"]]

      # generate synthetic indexes
      synthetic_treatment_idxs <- replicate(n = B, expr = sample.int(n = n_cells_curr_de, size = n_cells_curr_grna_group))

      # obtain vectors to pass to permutation test
      curr_expressions <- expressions[subset_vect]
      curr_fitted_means <- fitted_means[subset_vect]
      ground_truth_treatment_idxs <- seq(1, n_cells_curr_grna_group)

      # call the low-level association test function
      perm_out <- run_permutation_test(expressions = curr_expressions,
                                       fitted_means = curr_fitted_means,
                                       ground_truth_treatment_idxs = ground_truth_treatment_idxs,
                                       synthetic_treatment_idxs = synthetic_treatment_idxs,
                                       gene_theta = gene_theta,
                                       side = side,
                                       full_output = full_output)
    }) |> t() |> data.table::as.data.table() |>
      dplyr::mutate(grna_group = grna_groups, gene_id = gene_id)
  }) |> data.table::rbindlist()
}
