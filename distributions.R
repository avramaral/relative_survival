### Log Normal (LN)

# Hazard Function LN
hazLN <- function (N, time, mu, sigma, log = T, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, dlnorm(time[i], mu, sigma, log = T) - plnorm(time[i], mu, sigma, lower.tail = F, log.p = T))
  }
  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

# Cumulative Hazard Function LN
cumHazLN <- function (N, time, mu, sigma, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, - plnorm(time[i], mu, sigma, lower.tail = F, log.p = T))
  }
  res
}

### Log Logistic (LG)

# Hazard Function LG
hazLN <- function (N, time, mu, sigma, log = T, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, dlogis(log(time[i]), mu, sigma, log = T) - log(time[i]) - plogis(log(time[i]), mu, sigma, lower.tail = F, log.p = T))
  }
  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

# Cumulative Hazard Function LG
cumHazLN <- function (N, time, mu, sigma, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, - plogis(log(time[i]), mu, sigma, lower.tail = F, log.p = T))
  }
  res
}

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
