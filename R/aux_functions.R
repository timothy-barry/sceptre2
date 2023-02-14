#' Get gRNA index info
#'
#' Returns relevant information about the positions of the gRNAs among the cells.
#'
#' @param grna_group_assignments a vector of cell-to-gRNA assignments
#'
#' @return a list containing (i) the indexes of the cells to which each gRNA belongs, and (ii) the number of cells in which each gRNA belongs
get_grna_group_info <- function(grna_group_assignments, input_grna_groups) {
  unique_grnas <- c(input_grna_groups |> stats::setNames(input_grna_groups),
                    "non-targeting" = "non-targeting")
  # get the indices of each gRNA
  grna_specific_idxs <- lapply(unique_grnas, function(unique_grna) {
    which(grna_group_assignments == unique_grna)
   })
  # get the number of cells per gRNA; also record the number of NT cells
  n_cells_per_grna <- table(grna_group_assignments)[unique_grnas]
  return(list(grna_specific_idxs = grna_specific_idxs,
              n_cells_per_grna = n_cells_per_grna))
}


#' Get gRNA permutation indices
#'
#' Obtains the matrix of permutation indices for a given gRNA and choice of B
#'
#' @param n_cells_per_grna the "n_cells_per_grna" table of the grna_group_info list
#' @param unique_grna a gRNA
#' @param B the number of resamples to draw
#'
#' @return a matrix of permutation indices
get_grna_permutation_idxs <- function(n_cells_per_grna, unique_grna, B) {
  n_nt_cells <- n_cells_per_grna[["non-targeting"]]
  n_cells_curr_grna_group <- n_cells_per_grna[[unique_grna]]
  n_cells_curr_de <- n_cells_curr_grna_group + n_nt_cells
  synthetic_treatment_idxs <- replicate(n = B, expr = sample.int(n = n_cells_curr_de,
                                                                 size = n_cells_curr_grna_group))
}


#' Compute empirical p-value
#'
#' Computes an empirical permutation test p-value.
#'
#' @param z_star ground truth test statistic
#' @param z_null the null test statistics
#' @param side side of the tests
#'
#' @return the empirical p-value
compute_empirical_p_value <- function(z_star, z_null, side) {
  out_p <- switch(side,
                  "left" = mean(c(-Inf, z_null) <= z_star),
                  "right" = mean(c(Inf, z_null) > z_star),
                  "both" = 2 * min(mean(c(-Inf, z_null) <= z_star),
                                   mean(c(Inf, z_null) > z_star)))
  return(out_p)
}


#' Compute empirical p-value from gRNA-wise resul
#'
#' @param grna_wise_result data table giving the number of test statistics falling to the left of the original test statistic across batches (with genes in columns and rows in batches)
#' @param B total number of resamples
#' @param side sidedness of the test
#'
#' @return a vector of gene-wise p-values
compute_empirical_p_value_from_batch_result <- function(grna_wise_result, B, side) {
  if (is(grna_wise_result, "integer")) {
    left_tailed_p <- (grna_wise_result + 1)/B
  } else {
    left_tailed_p <- (apply(grna_wise_result, 2, sum) + 1)/B
  }
  switch(side,
         "left" = left_tailed_p,
         "right" = 1 - left_tailed_p,
         "both" = 2 * pmin(left_tailed_p, 1 - left_tailed_p))
}


#' Get target assignments via max operation
#'
#' @param grna_odm a grna ODM
#'
#' @return a vector of "targets" obtained by applying a column-wise max operation to the grna ODM;
#' the target is present as a column of the feature covariate matrix of the grna ODM.
#' @export
#'
#' @examples
#' \dontrun{
#' grna_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' get_target_assignments_via_max_op(grna_odm)
#' }
get_target_assignments_via_max_op <- function(grna_odm) {
  grna_feature_covariates <- grna_odm |> ondisc::get_feature_covariates()
  grna_feature_covariates$target[is.na(grna_feature_covariates$target)] <- "candidate"
  grna_to_target_map <- stats::setNames(row.names(grna_feature_covariates), grna_feature_covariates$target)
  grna_assignments <- get_grna_assignments_via_max_op(grna_odm)
  grna_targets <- names(grna_to_target_map)[match(x = grna_assignments, table = grna_to_target_map)]
  return(grna_targets)
}


#' Get grna assignments via max operation
#'
#' @param grna_odm a grna ODM;
#'
#' @return a vector of grna IDs obtained by applying a column-wise max operation to the grna ODM.
get_grna_assignments_via_max_op <- function(grna_odm) {
  ret <- grna_odm |>
    load_whole_odm() |>
    apply(MARGIN = 2, FUN = function(col) names(which.max(col))) |>
    unname()
  return(ret)
}


#' Load whole odm
#'
#' Loads data from disk into memory by cell.
#'
#' @param odm an ondisc_matrix object
#' @param csc_format load in cell-accessible, CSC format (TRUE) or feature-accessible, CSR format (FALSE)?
#'
#' @return an in-memory matrix
#' @export
load_whole_odm <- function(odm, csc_format = TRUE) {
  x <- odm@ondisc_matrix
  x_dim <- dim(x)
  index_on_cell <- csc_format
  subset_vector <- ondisc:::get_subset_vector(x, index_on_cell)
  if (identical(subset_vector, NA_integer_)) {
    subset_vector <- seq(1, if (index_on_cell) x_dim[2] else x_dim[1])
  }
  out <- ondisc:::return_spMatrix_from_index(x@h5_file, subset_vector,
                                             index_on_cell, x@logical_mat, x@underlying_dimension)
  second_subset <- ondisc:::get_subset_vector(x, !index_on_cell)
  if (!identical(second_subset, NA_integer_)) {
    out <- if (index_on_cell)
      out[second_subset, , drop = FALSE]
    else out[, second_subset, drop = FALSE]
  }
  row.names(out) <- ondisc::get_feature_ids(odm)
  colnames(out) <- ondisc::get_cell_barcodes(odm)
  return(out)
}
