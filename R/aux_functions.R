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


plot_fitted_density_result_row <- function(row, legend = TRUE) {
  z_null <- as.numeric(row[,grepl(pattern = "z_null_*", x = names(row))])
  z_star <- as.numeric(row[,"z_value"])
  dp <- row |>
    dplyr::select_if(names(row) %in% c("xi", "omega", "alpha", "nu"))  |>
    as.numeric()
  distribution <- if (length(dp) == 4) "ST" else "SN"
  p_out <- plot_fitted_density(dp = dp, z_null = z_null,
                               distribution = distribution,
                               z_star = z_star,
                               legend = legend)
  return(p_out)
}
