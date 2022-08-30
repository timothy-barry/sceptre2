#' Fit a skew distribution
#'
#' Fit a skew-t or skew-normal distribution on a vector of data.
#'
#' @param y a numeric vector of data.
#' @param max_it maximum number of iterations.
#'
#' @return a list containing (i) dp (the fitted model parameters), (ii) convergence success (TRUE or FALSE), and (iii) the number of iterations
#'
#' @examples
#' \dontrun{
#' y <- sn::rsn(n = 2500, dp = c(-.5, 1.5, 2))
#' fit_st <- fit_skew_t(y)
#' fit_sn <- fit_skew_normal(y)
#' plot_fitted_density(fit_st$dp, y, "ST")
#' plot_fitted_density(fit_sn$dp, y, "SN")
#'
#' y <- rpois(2000, 0.01)
#' fit_st <- fit_skew_t(y)
#' fit_sn <- fit_skew_normal(y)
#' }
fit_skew_t <- function(y, max_it = 70) {
  w <- rep(1, length(y))
  x <- matrix(1, nrow = length(y))
  fit <- sn::st.mple(x = x, y = y, w = w, control = list(iter.max = max_it))
  fit_info <- fit$opt.method
  return(list(dp = fit$dp, convergence = (fit_info$convergence == 0), n_iterations = fit_info$iterations))
}

#' @rdname fit_skew_t
fit_skew_normal <- function(y, max_it = 70) {
  w <- rep(1, length(y))
  x <- matrix(1, nrow = length(y))
  fit <- sn::sn.mple(x = x, y = y, w = w, control = list(iter.max = max_it))
  coef <- sn::cp2dp(cp = fit$cp, family = "SN")
  fit_info <- fit$opt.method
  return(list(dp = coef, convergence = (fit_info$convergence == 0), n_iterations = fit_info$iterations))
}


#' Compute skew distribution p-value
#'
#' Compute a p-value given a set of model parameters and a ground truth test statistic. Both skew-t and skew-normal options are available.
#'
#' @param dp skew-t or skew-normal distribution parameters
#' @param z_star ground truth test statistic
#' @param side sidedness of the test, one of "left", "right", and "both"
#'
#' @return a p-value
compute_skew_t_p_value <- function(dp, z_star, side) {
  switch(side,
         "left" = pmax(.Machine$double.eps, sn::pst(x = z_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
         "right" = pmax(.Machine$double.eps, 1 - sn::pst(x = z_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
         "both" = pmax(.Machine$double.eps, sn::pst(x = -abs(z_star), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
                         (1 - sn::pst(x = abs(z_star), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
         )
}


#' @rdname compute_skew_t_p_value
compute_skew_normal_p_value <- function(dp, z_star, side) {
  switch(side,
         "left" = pmax(.Machine$double.eps, sn::psn(x = z_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
         "right" = pmax(.Machine$double.eps, 1 - sn::psn(x = z_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
         "both" = pmax(.Machine$double.eps, sn::psn(x = -abs(z_star), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
                         (1 - sn::psn(x = abs(z_star), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
         )
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
                  "right" = mean(c(Inf, z_null) >= z_star),
                  "both" = 2 * min(mean(c(-Inf, z_null) <= z_star),
                                   mean(c(Inf, z_null) >= z_star)))
  return(out_p)
}


#' Compute KS test
#'
#' Runs a KS test on the resampled test statistics
#'
#' @param z_null the vector of null test statistics
#' @param dp the fitted parameters of the distribution
#' @param distribution the distribution, either "SN" (for skew-normal) or "ST" (for skew-t)
#'
#' @return a list containing (i) the KS test statistic, and (ii) the KS p-value
compute_ks_test <- function(z_null, dp, distribution) {
  if (distribution == "SN") {
    cdf <- function(x) sn::psn(x = x, dp = dp)
  } else if (distribution == "ST") {
    cdf <- function(x) sn::pst(x = x, dp = dp)
  }
  ret <- stats::ks.test(x = z_null, y = cdf)
  out <- stats::setNames(c(ret$statistic, ks_p_value = ret$p.value),
                         c("ks_stat", "ks_p_value"))
  return(out)
}


#' Plot fitted density
#'
#' Plots a fitted skew-normal or skew-t density superimposed over the histogram of resampled test statistics to which the density was fitted.
#'
#' @param dp vector of fitted skew-normal or skew-t parameters
#' @param z_null vector of resampled test statistics to which the model was fitted
#' @param distribution either "ST" (for skew-t) or "SN" (for skew-normal)
#' @param z_star (optional) ground truth test statistic
#'
#' @return a ggplot object
plot_fitted_density <- function(dp, z_null, distribution, z_star = NULL) {
  interval <- range(c(z_null, z_star)) + c(-0.25, 0.25)
  z <- seq(interval[1], interval[2], length.out = 1000)
  if (distribution == "ST") {
    fitted <- sn::dst(x = z, dp = dp)
    lab <- "Skew-t"
  } else if (distribution == "SN") {
    fitted <- sn::dsn(x = z, dp = dp)
    lab <- "Skew-normal"
  } else {
    stop("Distribution not recognized.")
  }
  df_curves <- data.frame(z = z, fitted = fitted, gaussian = stats::dnorm(z)) |>
    tidyr::gather("curve", "y", fitted, gaussian) |>
    dplyr::mutate(curve = factor(curve, levels = c("fitted","gaussian"), labels = c(lab, "N(0,1)")))
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = z, y = ..density..),
                            data = data.frame(z = z_null),
                            boundary = 0, colour = "black", fill = "lightgray", bins = 15) +
    ggplot2::geom_line(ggplot2::aes(x = z, y = y, group = curve, colour = curve, linetype = curve), data = df_curves) +
    ggplot2::scale_colour_manual(values = c("darkorchid2", "black"), name = "Null distribution") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = "Null distribution") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw() + ggplot2::xlab("") +
    (if (is.null(z_star)) NULL else ggplot2::geom_vline(xintercept = z_star, col = "darkred", lwd = 0.7)) +
    ggplot2::theme(legend.position = c(0.85, 0.8),
                   legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  return(p)
}
