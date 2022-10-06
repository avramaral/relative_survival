result_processing <- function (fit, model, time, X_tilde = matrix(data = NA, nrow = length(time), ncol = 0), X = matrix(data = NA, nrow = length(time), ncol = 0), verbose = T, ...) {
  
  if (!class(fit) == "stanfit") {
    stop("'fit' must be a 'stanfit' object.")
  }
  
  fitted_data <- extract(fit) 
  fix_coeff <- substring(text = model, first = c(4, 5), last = c(4, 5))
  
  if (fix_coeff[1] == fix_coeff[2]) {
    m_tilde <- check_n_cov(var = fitted_data$beta)
    m <- check_n_cov(var = fitted_data$beta)
  } else {
    m_tilde <- check_n_cov(var = fitted_data$alpha)
    m <- check_n_cov(var = fitted_data$beta)
  }

  if ((ncol(X_tilde) != m_tilde) | (ncol(X) != m)) {
    stop("Provide 'X' and 'X_tilde' with correct number of columns.")
  }
  
  name  <- substring(text = model, first = c(1, 4), last = c(3, 7))
  dist  <- gsub(pattern = "_", replacement = "", x = name[1])
  model <- name[2]
  
  spatial <- ifelse(test = as.logical(sum(str_detect(string = model, pattern = c("C", "D", "S", "T", "Y", "Z")))), yes = T, no = F)
  
  if (!spatial) {
    N_reg <- 1
  } else {
    if (model %in% c("ABCD", "ABXD", "ABCC", "ABDD", "XBXD", "AACC", "AADD", "BBCC", "BBDD")) {
      N_reg <- ncol(fitted_data$v)
    } else {
      N_reg <- ncol(fitted_data$u)
    }
  }
  
  N_samples <- length(fitted_data$lp__)
  
  excHaz    <- array(data = 0, dim = c(length(time), N_samples, N_reg))
  excCumHaz <- array(data = 0, dim = c(length(time), N_samples, N_reg))
  netSur    <- array(data = 0, dim = c(length(time), N_samples, N_reg))
  
  for (j in 1:N_reg) {
    
    if (verbose) { 
      cat(ifelse(test = spatial, yes = paste("\nRegion ", sprintf('%02d', j), "\n", sep = ""), no = "ALL\n"))
      progressbar <- txtProgressBar(min = 1, max = N_samples, initial = 1) 
    }
    
    # Vectorize this operation
    for (i in 1:N_samples) {
      if (fix_coeff[1] == fix_coeff[2]) {
        lp_tilde <- compute_lp(m = m_tilde, X = X_tilde, coeff = fitted_data$beta[i, ])
        lp       <- compute_lp(m = m, X = X, coeff = fitted_data$beta[i, ])
      } else {
        lp_tilde <- compute_lp(m = m_tilde, X = X_tilde, coeff = fitted_data$alpha[i, ])
        lp       <- compute_lp(m = m, X = X, coeff = fitted_data$beta[i, ])
      }
      
      part     <- add_re(fitted_data = fitted_data, model = model, lp_tilde = lp_tilde, lp = lp, i = i, j = j)
      lp_tilde <- part$lp_re_tilde
      lp       <- part$lp_re
      
      if (dist == "PGW") {
        excHaz[, i, j]    <- hazPGW(N = length(time), time = time * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i], log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazPGW(N = length(time), time = time * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i]) * exp(lp - lp_tilde)
      } else if (dist == "LN") {
        excHaz[, i, j]    <- hazLN(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i], log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazLN(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i]) * exp(lp - lp_tilde)
      } else if (dist == "LL") {
        excHaz[, i, j]    <- hazLL(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i], log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazLL(N = length(time), time = time * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i]) * exp(lp - lp_tilde)
      } else if (dist == "GG") {
        excHaz[, i, j]    <- hazGG(N = length(time), time = time * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i], log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazGG(N = length(time), time = time * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i]) * exp(lp - lp_tilde)
      } else if (dist == "GAM") {
        excHaz[, i, j]    <- hazGAM(N = length(time), time = time * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], log = F) * exp(lp)
        excCumHaz[, i, j] <- cumHazGAM(N = length(time), time = time * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i]) * exp(lp - lp_tilde)
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
