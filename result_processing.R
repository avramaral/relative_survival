result_processing <- function (model, fitted_data, N_samples, N_reg, distribution, time, X_tilde, X, verbose = T, ...) {
 
  m_tilde <- check_n_cov(var = fitted_data$alpha)
  m <- check_n_cov(var = fitted_data$beta)
  
  if ((ncol(X_tilde) != m_tilde) | (ncol(X) != m)) {
    stop("Provide the correct values for covariates.")
  }
  
  excHaz <- array(data = 0, dim = c(length(time), N_samples, N_reg))
  excCumHaz <- array(data = 0, dim = c(length(time), N_samples, N_reg))
  netSur <- array(data = 0, dim = c(length(time), N_samples, N_reg))
  
  for (j in 1:N_reg) {
    if (verbose) { 
      cat(paste("\nRegion ", sprintf('%02d', j), "\n", sep = ""))
      progressbar <- txtProgressBar(min = 1, max = N_samples, initial = 1) 
    }
    for (i in 1:N_samples) {
      lp_tilde <- compute_lp(m = m_tilde, X = X_tilde, coeff = fitted_data$alpha[i, ])
      lp <- compute_lp(m = m, X = X, coeff = fitted_data$beta[i, ])
      
      part <- add_re(model = model, lp_tilde = lp_tilde, lp = lp, random_effects = c(fitted_data$u_tilde[i, j], fitted_data$u[i, j]))
      lp_tilde <- part$lp_re_tilde
      lp <- part$lp_re
      
      if (distribution == "LN") {
        excHaz[, i, j] <- hazLN(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = exp(fitted_data$log_sigma[i]), log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazLN(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = exp(fitted_data$log_sigma[i])) * exp(lp - lp_tilde)
      } else if (distribution == "LL") {
        excHaz[, i, j] <- hazLL(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = exp(fitted_data$log_sigma[i]), log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazLL(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = exp(fitted_data$log_sigma[i])) * exp(lp - lp_tilde)
      } else if (distribution == "PGW") {
        excHaz[, i, j] <- hazPGW(N = length(time), time = time * exp(lp_tilde), eta = exp(fitted_data$log_eta[i]), nu = exp(fitted_data$log_nu[i]), theta = exp(fitted_data$log_theta[i]), log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazPGW(N = length(time), time = time * exp(lp_tilde), eta = exp(fitted_data$log_eta[i]), nu = exp(fitted_data$log_nu[i]), theta = exp(fitted_data$log_theta[i])) * exp(lp - lp_tilde)
      } else {
        stop("Choose a valid distribution.")
      }
      netSur[, i, j] <- exp(- excCumHaz[, i, j])
      if (verbose) { setTxtProgressBar(progressbar, i) }
    }
    if (verbose) { close(progressbar) } 
  }
  
  list(excHaz = excHaz, excCumHaz = excCumHaz, netSur = netSur)
}
