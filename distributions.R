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

### Log Logistic (LL)

# Hazard Function LG
hazLL <- function (N, time, mu, sigma, log = T, ...) {
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

# Cumulative Hazard Function LL
cumHazLL <- function (N, time, mu, sigma, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, - plogis(log(time[i]), mu, sigma, lower.tail = F, log.p = T))
  }
  res
}

### Generalized Gamma (GG)

# Hazard Function GG
hazGG <- function (N, time, eta, nu, theta, log = T, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, log(theta) - nu * log(eta) - lgamma(nu / theta) + (nu - 1) * log(time[i]) - (time[i] / eta) ^ theta - pgamma(q = time[i] ^ theta, shape = nu / theta, scale = eta ^ theta, log.p = T, lower.tail = F))
  }
  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

# Cumulative Hazard Function GG
cumHazGG <- function (N, time, eta, nu, theta, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, - pgamma(time[i] ^ theta, shape = nu / theta, scale = eta ^ theta, log = T, lower.tail = F))
  }
  res
}

### Gamma (GAM)

# Hazard Function GAM
hazGAM <- function (N, time, eta, nu, log = T, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, dgamma(t, shape = nu, scale = eta, log = T) - pgamma(t, shape = nu, scale = eta, lower.tail = F, log.p = T))
  }
  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

# Cumulative Hazard Function GAM
cumHazGAM <- function (N, time, eta, nu, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, - pgamma(t, shape = nu, scale = eta, lower.tail = F, log.p = T))
  }
  res
}
