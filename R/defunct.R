
#opt <- optim(dp, fn = st.pdev, gr = st.pdev.gh, method = opt.method,
#             control = control, x = x, y = y, w = w, fixed.nu = fixed.nu,
#             symmetr = symmetr, penalty = penalty.fn, trace = trace)
#my_st_obj <- function(dp, y) {
#  logL <- sum(dst(y, dp[1], dp[2], dp[3], dp[4], log = TRUE))
#  return(-2 * logL)
#}
#my_st_grad <- function(dp, y) {
#  xi <- dp[1]
#  omega <- dp[2]
#  alpha <- dp[3]
#  nu <- dp[4]
#  score <- numeric(4)
#  z <- (y - xi)/omega

#  nuz2 <- (nu + z^2)
#  loro.tau <- sqrt((nu + 1)/nuz2)
#  zt <- z * loro.tau
#  log.pdf <- dt(alpha * zt, nu + 1, log = TRUE)
#  log.cdf <- pt(alpha * zt, nu + 1, log.p = TRUE)
#  cdf <- exp(log.cdf)
#  loro.w <- exp(log.pdf - log.cdf)
#  tw <- loro.tau * loro.w
#  zwz2 <- z * (z^2 - 1) * loro.w/loro.tau
#  wi.beta <- z * loro.tau^2 - nu * alpha * tw/(nu + z^2)
#  score[1] <- sum(wi.beta)/omega
#  score[2] <- sum(-1 + zt^2 - alpha * nu * z * tw/(nu + z^2))/omega
#  score[3] <- sum(w * z * tw)
#  logTwz <- function(nu, alpha, z) {
#    r <- sqrt((nu + 1)/(nu + z^2))
#    pt(alpha * z * r, df = nu + 1, log.p = TRUE)
#  }
#  DlogTwz <- numDeriv::jacobian(logTwz, nu, z = z, alpha = alpha)
#  score[4] <- 0.5 * sum(-1/nu + digamma((nu + 1)/2) - digamma(nu/2) - log(1 + z^2/nu) + (nu + 1) * z^2/(nu * (nu + z^2)) + 2 * DlogTwz)
#  gradient <- (-2) * score
#  return(gradient)
#}
