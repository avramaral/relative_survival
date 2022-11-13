rneigh_matrix <- function (N_regions, lim_nodes = 10, silent = T, ...) {
  while (TRUE) {
    try({
      # Exponential degree distribution
      degrees <- sample(x = seq(from = 1, to = min(N_regions, lim_nodes)), size = N_regions, replace = T, prob = exp(-0.5 * seq(from = 1, to = min(N_regions, lim_nodes))))
      if (sum(degrees) %/% 2 != 0) { 
        degrees[1] <- degrees[1] + 1 
      }
      m <- as.matrix(as_adjacency_matrix(sample_degseq(out.deg = degrees, method = "vl")))
      break
    }, silent = silent)
  }
  m
}

# REF: https://github.com/atredennick/sageAbundance/
rICAR <- function (Q, ...) {
  e <- eigen(x = Q, symmetric = TRUE)
  v <- sqrt(ifelse(test = e$values > sqrt(.Machine$double.eps), yes = (1 / e$values), no = 0))
  sim <- e$vectors %*% diag(v) %*% rnorm(dim(Q)[1], 0, 1)
  X <- rep(x = 1, times = length(sim))
  if (sum(v == 0) == 2) {
    X <- cbind(X, seq(from = 1, to = length(sim)))
  }
  sim <- sim - X %*% solve(crossprod(x = X), crossprod(x = X, y = sim))
  sim
}

rIID <- function (sigma, N_regions, ...) {
  rnorm(n = N_regions, mean = 0, sd = sigma)
}

compute_re <- function (re, regions, ...) {
  r <- c()
  for (i in 1:length(regions)) {
    r <- c(r, re[regions[i]])
  }
  r
}

simulate_RS_MEGH <- function (dist, pars, X_tilde, X, alphas, betas, dep, re_tilde, re, idxs, ...) {
  
  compute_quantile <- function (dist, p, pars, ...) {
    if (dist == "PGW") {
      quant <- qPGW(p = p, eta = pars$eta, nu = pars$nu, theta = pars$theta)
    } else if (dist == "LN") {
      quant <- qLN(p = p, mu = pars$mu, sigma = pars$sigma)
    } else if (dist == "LL") {
      quant <- qLL(p = p, mu = pars$mu, sigma = pars$sigma)
    } else if (dist == "GG") {
      quant <- qGG(p = p, eta = pars$eta, nu = pars$nu, theta = pars$theta)
    } else if (dist == "GAM") {
      quant <- qGAM(p = p, eta = pars$eta, nu = pars$nu)
    } else {
      stop("Choose a valid distribution")
    }
    quant
  }
  
  u <- runif(n = nrow(X))
  x_tilde_alpha <- matrix(data = NA, nrow = nrow(X), ncol = 1)
  x_beta  <- matrix(data = NA, nrow = nrow(X), ncol = 1)
  for (i in 1:nrow(X)) {
    alpha <- alphas[[dep[i]]] 
    beta  <- betas[[dep[i]]]
    x_tilde_alpha[i, ] <- X_tilde[i, ] %*% alpha
    x_beta[i, ] <- X[i, ] %*% beta
  }
  p <- 1 - exp(log(1 - u) / exp(((x_beta) + compute_re(re = re, regions = desMat$gor)) - (x_tilde_alpha + compute_re(re = re_tilde, regions = desMat$gor))))
  
  quant <- list()
  for (i in 1:length(idxs)) { quant[[i]] <- compute_quantile(dist = dist, p = as.matrix(p[idxs[[i]], ]), pars = pars) }
  
  tmp <- matrix(data = 0, nrow = nrow(p), ncol = 1)
  for (j in 1:nrow(p)) {
    l <- which(sapply(idxs, function (e) is.element(j, e)))
    tmp[j, 1] <- quant[[l]][which(idxs[[l]] == j)]
  }
  quant <- tmp
  
  den <- exp(x_tilde_alpha + compute_re(re = re_tilde, regions = desMat$gor))
  
  (quant / den)
}

simulate_re <- function (struc = "ICAR", precision_tilde = 1, precision = 1, England = T, N_regions = 10, ...) {
  if (England) {
    W <- rbind(c(0, 1, 1, 0, 0, 0, 0, 0, 0),
               c(1, 0, 1, 1, 1, 0, 0, 0, 0),
               c(1, 1, 0, 1, 0, 0, 0, 0, 0),
               c(0, 1, 1, 0, 1, 1, 0, 1, 0),
               c(0, 1, 0, 1, 0, 0, 0, 1, 1),
               c(0, 0, 0, 1, 0, 0, 1, 1, 0),
               c(0, 0, 0, 0, 0, 1, 0, 1, 0),
               c(0, 0, 0, 1, 1, 1, 1, 0, 1),
               c(0, 0, 0, 0, 1, 0, 0, 1, 0)) 
  } else {
    W <- rneigh_matrix(N_regions = N_regions)
  }
  D <- diag(rowSums(W))
  if (struc == "ICAR") {
    re_tilde <- rICAR(Q = precision_tilde * (D - W))
    re <- rICAR(Q = precision * (D - W))  
  } else if (struc == "IID") {
    re_tilde <- as.matrix(rIID(sigma = sqrt(1 / precision), N_regions = nrow(W)))
    re <- as.matrix(rIID(sigma = sqrt(1 / precision), N_regions = nrow(W)))
  } else if (struc == "NONE") {
    re_tilde <- as.matrix(rep(x = 0, times = nrow(W)))
    re <- as.matrix(rep(x = 0, times = nrow(W)))
  }
  list(re_tilde = re_tilde, re = re)
}
