#' Compute Fisher's exact p-value
#'
#' @param treatment_expressions the vector of expressions for the treatment cells (i.e., the cells receiving a targeting gRNA)
#' @param control_expressions the vector of expressions for the control cells (i.e., the cells receiving an NT)
#'
#' @return p-value from the Fisher exact test for two-way contingency tables
compute_fisher_exact_p_value <- function(n_treatment_cells_with_expression,
                                         n_treatment_cells_without_expression,
                                         n_control_cells_with_expression,
                                         n_control_cells_without_expression) {
  contingency_table <- matrix(c(n_treatment_cells_with_expression,
                                n_treatment_cells_without_expression,
                                n_control_cells_with_expression,
                                n_control_cells_without_expression),
                              nrow = 2, byrow = TRUE)
  fit <- stats::fisher.test(contingency_table)
  return(fit$p.value)
}


#' Get contingency table
#'
#' @param curr_expressions the vector of current expressions
#' @param ground_truth_treatment_idxs index vector indicating which of the cells are the treatment cells
#'
#' @return returns a vector containing (i) the number of "treatment" cells, (ii) the number of "control" cells, (iii) the number of treatment cells with nonzero expressions, and (iv) the number of control cells with nonzero expressions.
get_contingency_table <- function(curr_expressions, ground_truth_treatment_idxs) {
  treatment_expressions <- curr_expressions[ground_truth_treatment_idxs]
  control_expressions <- curr_expressions[-ground_truth_treatment_idxs]

  n_treatment_cells <- length(treatment_expressions)
  n_treatment_cells_with_expression <- sum(treatment_expressions >= 1)
  n_treatment_cells_without_expression <- n_treatment_cells - n_treatment_cells_with_expression

  n_control_cells <- length(control_expressions)
  n_control_cells_with_expression <- sum(control_expressions >= 1)
  n_control_cells_without_expression <- n_control_cells - n_control_cells_with_expression

  return(c(n_treatment_cells_with_expression = n_treatment_cells_with_expression,
         n_treatment_cells_without_expression = n_treatment_cells_without_expression,
         n_control_cells_with_expression = n_control_cells_with_expression,
         n_control_cells_without_expression = n_control_cells_without_expression))
}
