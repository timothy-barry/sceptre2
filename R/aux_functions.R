#' Get gRNA index info
#'
#' Returns relevant information about the positions of the gRNAs among the cells.
#'
#' @param grna_group_assignments a vector of cell-to-gRNA assignments
#'
#' @return a list containing (i) the indexes of the cells to which each gRNA belongs, and (ii) the number of cells in which each gRNA belongs
get_grna_group_info <- function(grna_group_assignments, input_grna_groups, B) {
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


get_grna_permutation_idxs <- function(n_cells_per_grna, unique_grna, B) {
  set.seed(4)
  n_nt_cells <- n_cells_per_grna[["non-targeting"]]
  n_cells_curr_grna_group <- n_cells_per_grna[[unique_grna]]
  n_cells_curr_de <- n_cells_curr_grna_group + n_nt_cells
  synthetic_treatment_idxs <- replicate(n = B, expr = sample.int(n = n_cells_curr_de,
                                                                 size = n_cells_curr_grna_group))
}


plot_fitted_density_result_row <- function(row, n_bins = 15, legend = TRUE) {
  z_null <- as.numeric(row[,grepl(pattern = "z_null_*", x = names(row))])
  z_star <- as.numeric(row[,"z_value"])
  dp <- row |>
    dplyr::select_if(names(row) %in% c("xi", "omega", "alpha", "nu"))  |>
    as.numeric()
  distribution <- if (length(dp) == 4) "ST" else "SN"
  p_out <- plot_fitted_density(dp = dp, z_null = z_null,
                               distribution = distribution,
                               z_star = z_star,
                               legend = legend,
                               n_bins = n_bins)
  return(p_out)
}
