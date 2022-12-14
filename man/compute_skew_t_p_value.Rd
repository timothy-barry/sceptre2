% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skew_distributions.R
\name{compute_skew_t_p_value}
\alias{compute_skew_t_p_value}
\alias{compute_skew_normal_p_value}
\title{Compute skew distribution p-value}
\usage{
compute_skew_t_p_value(dp, z_star, side)

compute_skew_normal_p_value(dp, z_star, side)
}
\arguments{
\item{dp}{skew-t or skew-normal distribution parameters}

\item{z_star}{ground truth test statistic}

\item{side}{sidedness of the test, one of "left", "right", and "both"}
}
\value{
a p-value
}
\description{
Compute a p-value given a set of model parameters and a ground truth test statistic. Both skew-t and skew-normal options are available.
}
