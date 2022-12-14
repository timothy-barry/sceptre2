% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contingency_table.R
\name{get_contingency_table}
\alias{get_contingency_table}
\title{Get contingency table}
\usage{
get_contingency_table(curr_expressions, ground_truth_treatment_idxs)
}
\arguments{
\item{curr_expressions}{the vector of current expressions}

\item{ground_truth_treatment_idxs}{index vector indicating which of the cells are the treatment cells}
}
\value{
returns a vector containing (i) the number of "treatment" cells, (ii) the number of "control" cells, (iii) the number of treatment cells with nonzero expressions, and (iv) the number of control cells with nonzero expressions.
}
\description{
Get contingency table
}
