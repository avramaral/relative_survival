### Power Generalized Weibull (PGW)

# Hazard Function PGW
hazPGW <- function (N, time, eta, nu, theta, log = T, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, log(nu) - log(theta) - nu * log(eta) + (nu - 1) * log(time[i]) + ((1 / theta) - 1) * log(1 + (time[i] / eta) ^ nu))
  }
  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

# Cumulative Hazard Function PGW
cumHazPGW <- function (N, time, eta, nu, theta, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, - 1 + (1 + (time[i] / eta) ^ nu) ^ (1 / theta))
  }
  res
}