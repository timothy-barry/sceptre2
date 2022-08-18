#' Fit skew t
#'
#' Fits a skew-t distribution to a vector of data y.
#'
#' @param y a numeric vector of data.
#' @param max_it maximum number of iterations.
#'
#' @return a list containing (i) dp, the fitted parameters and (ii) convergence success (TRUE or FALSE)
#' @export
#'
#' @examples
#' y <- sn::rst(n = 2500, dp = c(0, 0.5, 2, 20))
#' fit <- fit_skew_t(y)
#' plot_skew_t(fit$dp, y)
#'
#' y <- rpois(2000, 0.01)
#' fit <- fit_skew_t(y)
#'
#' sn::selm(y ~ 1, family = "ST")@param$dp
fit_skew_t <- function(y, max_it = 70) {
  w <- rep(1, length(y))
  x <- matrix(1, nrow = length(y))
  init_guess <- as.numeric(sn:::mst.prelimFit(x = x, y = y, w = w, quick = TRUE)$dp)
  tiny <- (.Machine$double.eps)^(0.25)
  low.dp <- c(-Inf, tiny, -Inf, tiny)
  high.dp <- rep(Inf, 4)
  opt <- nlminb(init_guess, objective = sn:::st.pdev, gradient = sn:::st.pdev.gh,
                lower = low.dp, upper = high.dp, control = list(iter.max = max_it),
                x = x, y = y, w = w)
  return(list(dp = opt$par, convergence = (opt$convergence == 0), n_iterations = opt$iterations))
}


#' Plot skew-t
#'
#' @param dp vector of fitted skew-t parameters
#' @param z_null vector of resampled test statistics to which skew-t was fitted
#' @param z_star (optional) ground truth test statistic
#'
#' @return a ggplot object
plot_skew_t <- function(dp, z_null, z_star = NULL) {
  interval <- range(c(z_null, z_star)) + c(-0.25, 0.25)
  z <- seq(interval[1], interval[2], length.out = 1000)
  df_curves <- data.frame(z = z, fitted = sn::dst(x = z, dp = dp), gaussian = stats::dnorm(z)) |>
    tidyr::gather("curve", "y", fitted, gaussian) |>
    dplyr::mutate(curve = factor(curve, levels = c("fitted","gaussian"), labels = c("Skew-t","N(0,1)")))
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
