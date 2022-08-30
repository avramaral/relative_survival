hazPGW <- function (time, eta, nu, theta) {
  N <- length(time)
  result <- c()
  for (i in 1:N) { 
    result <- c(result, nu / (theta * eta ^ nu) * time[i] ^ (nu - 1) * (1 + (time[i] / eta) ^ nu) ^ ((1 / theta) - 1))
  }
  result
}

cumHazPGW <- function (time, eta, nu, theta) {
  N <- length(time)
  result <- c()
  for (i in 1:N) { 
    result <- c(result, - 1 + (1 + (time[i] / eta) ^ nu) ^ (1 / theta))
  }
  result
}

estSurv <- function (cumHaz) {
  exp(-cumHaz)
}

linear_predictor <- function (X, fixed_coeff, region, random_effect) {
  if (ncol(X) != length(fixed_coeff)) {
    stop("Dimension mismatch.")
  } else {
    lp <- X %*% fixed_coeff + random_effect[as.integer(region)]
  }
  lp
}

