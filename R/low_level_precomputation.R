#' Run response precomputation
#'
#' Runs precomputation on a response
#'
#' @param expressions the numeric vector of response expressions
#' @param covariate_matrix the covariate matrix on which to regress (NOTE: should contain an interecept term)
#'
#' @return a named vector of fitted coefficients, alongside the fitted size parameter (named "response_theta")
run_response_precomputation_low_level <- function(expressions, covariate_matrix) {
  # second backup: method of moments
  backup_2 <- function(pois_fit) {
    MASS::theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
  }

  # first backup: MLE on poisson reg
  backup <- function() {
    pois_fit <- stats::glm(expressions ~ . + 0, data = covariate_matrix, family = stats::poisson())
    response_theta <- tryCatch({
      MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
    }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
    response_theta <- max(response_theta, 0.1)
    fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(response_theta), data = covariate_matrix)
    fitted_coefs <- stats::coef(fit_nb)
    return(c(fitted_coefs = fitted_coefs, response_theta = response_theta))
  }

  # try to fit a negative binomial GLM with unknown dispersion
  result <- tryCatch({
    fit_nb_init <- MASS::glm.nb(formula = expressions ~ . + 0, data = covariate_matrix)
    response_theta <- max(fit_nb_init$theta, 0.1)
    fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(response_theta), data = covariate_matrix)
    fitted_coefs <- stats::coef(fit_nb)
    return(c(fitted_coefs, response_theta = response_theta))
  }, error = function(e) backup(), warning = function(w) backup())

  return(result)
}


#' Run grna precomputation
#'
#' Runs precomputation on a vector of grna indicators
#'
#' @param indicators the binary vector of grna indicators
#' @param covariate_matrix the covariate matrix on which to regress (NOTE: should contain an interecept term)
#'
#' @return a named vector of fitted coefficients
run_grna_precomputation_low_level <- function(indicators, covariate_matrix) {
  fit <- stats::glm(formula = indicators ~ . + 0, family = "binomial", data = covariate_matrix)
  ret <- stats::coef(fit)
  return(ret)
}


#' Get pieces from response precomp
#'
#' @param response_precomp a response precomp vector
#' @param global_cell_covariates a global cell covariates matrix
#'
#' @return a list containing (i) response_theta, and (ii) the fitted linear components
get_pieces_from_response_precomp <- function(response_precomp, global_cell_covariates) {
  response_fitted_coefs <- response_precomp[-length(response_precomp)]
  response_theta <- response_precomp[[length(response_precomp)]]
  fitted_means <- exp((as.matrix(global_cell_covariates) %*% response_fitted_coefs)[,1])
  return(list(response_theta = response_theta, fitted_means = fitted_means))
}
