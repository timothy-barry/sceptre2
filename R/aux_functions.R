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
  left_tailed_p <- (apply(grna_wise_result, 2, sum) + 1)/B
  switch(side,
         "left" = left_tailed_p,
         "right" = 1 - left_tailed_p,
         "both" = 2 * pmin(left_tailed_p, 1 - left_tailed_p))

}
