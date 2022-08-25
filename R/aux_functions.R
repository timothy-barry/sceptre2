#' Get gRNA index info
#'
#' Returns relevant information about the positions of the gRNAs among the cells.
#'
#' @param grna_group_assignments a vector of cell-to-gRNA assignments
#'
#' @return a list containing (i) the indexes of the cells to which each gRNA belongs, and (ii) the number of cells in which each gRNA belongs
get_grna_group_info <- function(grna_group_assignments) {
  unique_grnas <- unique(grna_group_assignments)
  unique_grnas <- stats::setNames(unique_grnas, unique_grnas)
  grna_specific_idxs <- lapply(unique_grnas, function(unique_grna) {
    which(grna_group_assignments == unique_grna)
   })
  n_cells_per_grna <- table(grna_group_assignments)
  return(list(grna_specific_idxs = grna_specific_idxs,
              n_cells_per_grna = n_cells_per_grna))
}


compute_empirical_p_value_result_row <- function(row, side = "both") {
  z_null <- row[grepl(pattern = "z_null_*", x = names(row))]
  z_star <- row["z_value"]
  compute_empirical_p_value(z_star, z_null, side)
}


plot_fitted_density_result_row <- function(row) {
  z_null <- row[grepl(pattern = "z_null_*", x = names(row))]
  z_star <- row["z_value"]
  dp <- row[c("xi", "omega", "alpha", "nu")]
  dp <- dp[!is.na(dp)]
  distribution <- if (length(dp) == 4) "ST" else "SN"
  p_out <- plot_fitted_density(dp = dp, z_null = z_null,
                               distribution = distribution, z_star = z_star)
  p_out
}
