data_stan <- function (data, model, cov_tilde = c(), cov = c(), intercept_tilde = F, intercept = F, pop.haz = "pop.haz", time = "time", obs = "cens", region = "region", adj_info = list(), ...) {
  
  observed <- which(data[[obs]] == 1)
  N <- nrow(data)
  N_obs <- length(observed)
  pre_computed <- list(observed = observed, N = N, N_obs = N_obs)
  
  if (model == 1) {
    data <- data_modelX(data = data, pre_computed = pre_computed, cov_tilde = cov_tilde, cov = cov, intercept_tilde = intercept_tilde, intercept = intercept, pop.haz = pop.haz, time = time, region = region, adj_info  = adj_info)
  } else if (model == 2) {
    data <- data_modelX(data = data, pre_computed = pre_computed, cov_tilde = cov_tilde, cov = cov, intercept_tilde = intercept_tilde, intercept = intercept, pop.haz = pop.haz, time = time, region = region, adj_info  = adj_info)
  } else if (model == 3) {
    data <- data_modelX(data = data, pre_computed = pre_computed, cov_tilde = cov_tilde, cov = cov, intercept_tilde = intercept_tilde, intercept = intercept, pop.haz = pop.haz, time = time, region = region, adj_info  = adj_info)
  } else if (model == 4) {
    data <- data_modelY(data = data, pre_computed = pre_computed,                        cov = cov,                                    intercept = intercept, pop.haz = pop.haz, time = time, region = region, adj_info  = adj_info)
  } else if (model == 5) {
    data <- data_modelY(data = data, pre_computed = pre_computed,                        cov = cov,                                    intercept = intercept, pop.haz = pop.haz, time = time, region = region, adj_info  = adj_info)
  } else if (model == 6) {
    data <- data_model6(data = data, pre_computed = pre_computed, cov_tilde = cov_tilde, cov = cov, intercept_tilde = intercept_tilde, intercept = intercept, pop.haz = pop.haz, time = time)
  } else if (model == 7) {
    data <- data_modelZ(data = data, pre_computed = pre_computed,                        cov = cov,                                    intercept = intercept, pop.haz = pop.haz, time = time)
  } else if (model == 8) {
    data <- data_modelZ(data = data, pre_computed = pre_computed,                        cov = cov,                                    intercept = intercept, pop.haz = pop.haz, time = time)
  } else if (model == 9) {
    data <- data_model9(data = data, pre_computed = pre_computed, cov_tilde = cov_tilde,            intercept_tilde = intercept_tilde,                        pop.haz = pop.haz, time = time)
  } else {
    stop("Select a valid model.")
  }
  data
}

design_matrix <- function (data, N, cov, intercept, ...) {
  X_names <- cov
  int <- data.frame()[1:N, ]
  if (intercept) { int <- rep(1, N) }
  X <- as.matrix(cbind(int, data[X_names]))
  rownames(X) <- NULL
  if (intercept) { colnames(X) <- c("int", X_names) } else { colnames(X) <- X_names }
  M <- ncol(X)
  list(X = X, M = M)
}

data_modelX <- function (data, pre_computed, cov_tilde, cov, intercept_tilde, intercept, pop.haz, time, region, adj_info, ...) { # Works for models 1, 2, and 3
  
  if (length(adj_info) == 0) {
    stop("Provide 'adj_info.'")
  }
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov_tilde, intercept = intercept_tilde)
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, intercept = intercept)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X, X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]))
}

data_modelY <- function (data, pre_computed, cov, intercept, pop.haz, time, region, adj_info, ...) {
  
  if (length(adj_info) == 0) {
    stop("Provide 'adj_info.'")
  }
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, intercept = intercept)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]))
}

data_model6 <- function (data, pre_computed, cov_tilde, cov, intercept_tilde, intercept, pop.haz, time, ...) {
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov_tilde, intercept = intercept_tilde)
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, intercept = intercept)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X, X = dm$X)
}

data_modelZ <- function (data, pre_computed, cov, intercept, pop.haz, time, ...) {
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, intercept = intercept)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X)
}

data_model9 <- function (data, pre_computed, cov_tilde, intercept_tilde, pop.haz, time, ...) {
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov_tilde, intercept = intercept_tilde)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X)
}
