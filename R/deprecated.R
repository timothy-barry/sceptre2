# the old precomputation function that used VGAM to fit the GLM.
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
      fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(response_theta), data = covariate_matrix)
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
    fit_nb_init <- MASS::glm.nb(formula = expressions ~ . + 0, data = covariate_matrix)
    response_theta <- min(max(fit_nb_init$theta, 0.1), 1000)
    fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(response_theta), data = covariate_matrix)
    fitted_coefs <- stats::coef(fit_nb)
    list(precomp_str = "nb:mass", precomp = c(fitted_coefs, response_theta = response_theta))
  }, error = function(e) backup(), warning = function(w) backup())

  return(result)
}
