result_processing <- function (fit, model, time, X_tilde, X, region = NULL, verbose = T, ...) {
  
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
  
  if (spatial) {
    if (model %in% c("ABCD", "ABXD", "ABCC", "ABDD", "XBXD", "AACC", "AADD", "BBCC", "BBDD")) {
      N_reg <- ncol(fitted_data$v)
    } else {
      N_reg <- ncol(fitted_data$u)
    }
  }
  
  N_sim <- length(fitted_data$lp__)
  
  excHaz    <- array(data = 0, dim = c(nrow(X), N_sim, length(time)))
  excCumHaz <- array(data = 0, dim = c(nrow(X), N_sim, length(time)))
  netSur    <- array(data = 0, dim = c(nrow(X), N_sim, length(time)))
  
  if (verbose) { progressbar <- txtProgressBar(min = 1, max = N_sim, initial = 1) }
  
  for (i in 1:N_sim) {
    
    if (fix_coeff[1] == fix_coeff[2]) {
      lp_tilde <- compute_lp(m = m_tilde, X = as.matrix(X_tilde), coeff = fitted_data$beta[i, ])
      lp       <- compute_lp(m = m, X = as.matrix(X), coeff = fitted_data$beta[i, ])
    } else {
      lp_tilde <- compute_lp(m = m_tilde, X = as.matrix(X_tilde), coeff = fitted_data$alpha[i, ])
      lp       <- compute_lp(m = m, X = as.matrix(X), coeff = fitted_data$beta[i, ])
    }
    
    part     <- add_re_mod(fitted_data = fitted_data, model = model, lp_tilde = lp_tilde, lp = lp, i = i, j = region)
    lp_tilde <- part$lp_re_tilde
    lp       <- part$lp_re
    
    compute_haz <- function (t, dist, ...) {
      if (dist == "PGW") {
        haz <- hazPGW(N = nrow(X), time = time[t] * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i], log = F) * exp(lp)
      } else if (dist == "LN") {
        haz <-  hazLN(N = nrow(X), time = time[t] * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i], log = F) * exp(lp)
      } else if (dist == "LL") {
        haz <-  hazLL(N = nrow(X), time = time[t] * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i], log = F) * exp(lp)
      } else if (dist == "GG") {
        haz <-  hazGG(N = nrow(X), time = time[t] * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i], log = F) * exp(lp)
      } else if (dist == "GAM") {
        haz <- hazGAM(N = nrow(X), time = time[t] * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], log = F) * exp(lp)
      }
      haz
    }
    
    compute_cumHaz <- function (t, dist, ...) {
      if (dist == "PGW") {
        cumHaz <- cumHazPGW(N = nrow(X), time = time[t] * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i]) * exp(lp - lp_tilde)
      } else if (dist == "LN") {
        cumHaz <-  cumHazLN(N = nrow(X), time = time[t] * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i]) * exp(lp - lp_tilde)
      } else if (dist == "LL") {
        cumHaz <-  cumHazLL(N = nrow(X), time = time[t] * exp(lp_tilde), mu = fitted_data$mu[i], sigma = fitted_data$sigma[i]) * exp(lp - lp_tilde)
      } else if (dist == "GG") {
        cumHaz <-  cumHazGG(N = nrow(X), time = time[t] * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i], theta = fitted_data$theta[i]) * exp(lp - lp_tilde)
      } else if (dist == "GAM") {
        cumHaz <- cumHazGAM(N = nrow(X), time = time[t] * exp(lp_tilde), eta = fitted_data$eta[i], nu = fitted_data$nu[i]) * exp(lp - lp_tilde)
      }
      cumHaz
    }
    
    if (dist == "PGW") {
      excHaz[, i, ]    <- sapply(X = 1:length(time), compute_haz, dist = "PGW")
      excCumHaz[, i, ] <- sapply(X = 1:length(time), compute_cumHaz, dist = "PGW")
    } else if (dist == "LN") {
      excHaz[, i, ]    <- sapply(X = 1:length(time), compute_haz, dist = "LN")
      excCumHaz[, i, ] <- sapply(X = 1:length(time), compute_cumHaz, dist = "LN")
    } else if (dist == "LL") {
      excHaz[, i, ]    <- sapply(X = 1:length(time), compute_haz, dist = "LL")
      excCumHaz[, i, ] <- sapply(X = 1:length(time), compute_cumHaz, dist = "LL")
    } else if (dist == "GG") {
      excHaz[, i, ]    <- sapply(X = 1:length(time), compute_haz, dist = "GG")
      excCumHaz[, i, ] <- sapply(X = 1:length(time), compute_cumHaz, dist = "GG")
    } else if (dist == "GAM") {
      excHaz[, i, ]    <- sapply(X = 1:length(time), compute_haz, dist = "GAM")
      excCumHaz[, i, ] <- sapply(X = 1:length(time), compute_cumHaz, dist = "GAM")
    } else {
      stop("Choose a valid distribution.")
    }
    
    netSur[, i, ] <- exp(- excCumHaz[, i, ])
    
    if (verbose) { setTxtProgressBar(progressbar, i) }
    
  }
  
  if (verbose) { close(progressbar) }
  
  list(excHaz = excHaz, excCumHaz = excCumHaz, netSur = netSur)
}

#################
# Grouping step #
#################

compute_summary <- function (x, ...) {
  x_mean <- apply(X = x, MARGIN = 3, FUN = mean)
  x_Q025 <- apply(X = x, MARGIN = 3, FUN = quantile, probs = c(0.025))
  x_Q975 <- apply(X = x, MARGIN = 3, FUN = quantile, probs = c(0.975))
  list(mean = x_mean, Q025 = x_Q025, Q975 = x_Q975)
}

group_marginal_quantities <- function (group, excHaz, excCumHaz, netSur, ...) {
  
  n_cat <- length(unique(group))
  
  excHaz_summary    <- list()
  excCumHaz_summary <- list()
  netSur_summary    <- list()
  
  progressbar <- txtProgressBar(min = 1, max = n_cat, initial = 1)
  
  for (i in 1:n_cat) {
    idx <- which(group == i)
    
    excHaz_summary[[i]]    <- compute_summary(x = excHaz[idx, , ])
    excCumHaz_summary[[i]] <- compute_summary(x = excCumHaz[idx, , ])
    netSur_summary[[i]]    <- compute_summary(x = netSur[idx, , ])
    
    setTxtProgressBar(progressbar, i)
  }
  close(progressbar)
  
  list(excHaz_summary = excHaz_summary, excCumHaz_summary = excCumHaz_summary, netSur_summary = netSur_summary)
}
