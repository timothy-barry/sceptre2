perform_association_test_lowmoi_odm <- function(mm_odm, grna_group_assignment, nt_cell, gene_grna_group_pairs, B, full_output, side) {
  # obtain gene odm
  gene_odm <- mm_odm |> ondisc::get_modality("gene")

  # get ordered gene-grna pairs, as well as the vector of unique gene ids
  gene_grna_group_pairs <- gene_grna_group_pairs |> dplyr::arrange(gene_id, grna_group)
  gene_ids <- as.character(unique(gene_grna_group_pairs$gene_id))

  # get the global cell covariates
  global_cell_covariates <- mm_odm |> get_cell_covariates()

  # tabulate the number of cells per gRNA
  grna_group_tab <- table(grna_group_assignment)
  n_nt_cells <- grna_group_tab[["non-targeting"]]
  n_test_cells_per_grna_group <- grna_group_tab[names(grna_group_tab) != "non-targeting"] + n_nt_cells

  # loop over the gene ids
  for (gene_id in gene_ids) {
    # load expression data
    expressions <- as.numeric(gene_odm[[gene_id,]])

    # perform the gene precomputation
    gene_precomp <- run_gene_precomputation_low_level(expressions = expressions[nt_cell],
                                                      covariate_matrix = global_cell_covariates[nt_cell,])

    # compute the fitted values of the regression
    pieces <- get_pieces_from_gene_precomp(gene_precomp = gene_precomp,
                                           global_cell_covariates = global_cell_covariates)
    gene_theta <- pieces$gene_theta
    fitted_means <- pieces$fitted_means

    # iterate over the gRNAs for a given gene, generating the permutations and fitting the model
    grna_groups <- dplyr::filter(gene_grna_group_pairs, gene_id == !!gene_id) |>
      dplyr::pull(grna_group) |>
      as.character()

    for (grna_group in grna_groups) {
      # for each grna_group, generate a matrix of permutation indexes
      synthetic_treatment_idxs <- replicate(n = B, expr = sample.int(n = n_test_cells_per_grna_group[[grna_group]],
                                                                     size = grna_group_tab[[grna_group]]))


      # subset the expressions vector to filter for "treatment" or "control" cells
      # NOTE: Consider factoring out this step from the loop
      curr_idxs <- grna_group_assignment %in% c("non-targeting", grna_group)
      subsetted_grna_group_assignments <- grna_group_assignments[curr_idxs]
      subsetted_expressions <- expressions[curr_idxs]
      subsetted_fitted_means <- fitted_means[curr_idxs]
      subsetted_ground_truth_treatment_idxs <- which(subsetted_grna_group_assignments == grna_group)

      # call the low-level association test function
      run_permutation_test(expressions = subsetted_expressions,
                           fitted_means = subsetted_fitted_means,
                           ground_truth_treatment_idxs = subsetted_ground_truth_treatment_idxs,
                           synthetic_treatment_idxs = synthetic_treatment_idxs,
                           gene_theta = gene_theta,
                           side = side,
                           full_output = full_output)
    }
  }
}
