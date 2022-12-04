#' Run response precomputation
#'
#' Runs precomputation on a response
#'
#' @param expressions the numeric vector of response expressions
#' @param covariate_matrix the covariate matrix on which to regress (NOTE: should contain an interecept term)
#'
#' @return a list containing the following elements:
#' (i) "precomp_str": a string summarizing the method for fitting the GLM and the size parameter
#' (ii) "fitted_coefs": a vector of fitted coefficients; the final entry of this vector is the fitted theta
run_response_precomputation_low_level <- function(expressions, covariate_matrix) {
  # backup: return fitted coefficients from Poisson regression
  backup_3 <- function(pois_fit, pois_warn) {
    list(fitted_coef_str = paste0("pois", if (pois_warn) "_(warn)" else NULL),
         fitted_coefs = stats::coef(pois_fit))
  }

  # Backup: method of moments
  backup_2 <- function(pois_fit) {
    list(theta_fit_str = "resid_mm",
         MASS::theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual))
  }

  #Backup: MLE on poisson reg
  backup <- function() {
    # fit Poisson model, tracking warning
    pois_warn <- FALSE
    wHandler <- function(w) {pois_warn <<- TRUE; invokeRestart("muffleWarning")}
    withCallingHandlers(expr = {
      pois_fit <- stats::glm(expressions ~ . + 0, data = covariate_matrix, family = stats::poisson())
    }, warning = wHandler)

    # get theta; save theta itself (in response_theta) and theta_fit_string (i.e., procedure used to fit theta)
    response_theta_list <- tryCatch({
      list(theta_fit_str = "resid_mle",
           response_theta = MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1])
    }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
    theta_fit_str <- response_theta_list$theta_fit_str
    response_theta <- response_theta_list$response_theta
    response_theta <- min(max(response_theta, 0.1), 1000)

    # obtain the fitted coefficients
    fitted_coefs_list <- tryCatch({
      fit_nb <- stats::glm(formula = expressions ~ . + 0, family = MASS::negative.binomial(response_theta), data = covariate_matrix)
      list(fitted_coef_str = "nb",
           fitted_coefs = stats::coef(fit_nb))
    }, error = function(e) backup_3(pois_fit, pois_warn), warning = function(w) backup_3(pois_fit, pois_warn))
    fitted_coefs <- fitted_coefs_list$fitted_coefs
    precomp <- c(fitted_coefs, response_theta = response_theta)
    fitted_coef_str <- fitted_coefs_list$fitted_coef_str

    # construct the precomp string
    precomp_str <- paste0(fitted_coef_str, ":", theta_fit_str)

    list(precomp_str = precomp_str,
         precomp = c(fitted_coefs, response_theta = response_theta))
  }

  # try to fit a negative binomial GLM with unknown dispersion
  result <- tryCatch({
    fit_nb <- MASS::glm.nb(formula = expressions ~ . + 0, data = covariate_matrix)
    response_theta <- min(max(fit_nb_init$theta, 0.1), 1000)
    fitted_coefs <- stats::coef(fit_nb)
    list(precomp_str = "nb:mass",
         precomp = c(fitted_coefs, response_theta = response_theta))
  }, error = function(e) backup(), warning = function(w) backup())

  return(result)
}


#' Run grna precomputation
#'
#' Runs precomputation on a vector of grna indicators.
#' This function is used for conditional resampling tests only.
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
  mu_hats <- exp(as.numeric((global_cell_covariates %*% response_fitted_coefs)))
  return(list(response_theta = response_theta,
              mu_hats = mu_hats))
}
