data_stan <- function(data, model, cov.tilde = c(), cov = c(), nonlinear = c(), df = 3, time = "time", obs = "obs", pop.haz = "pop.haz", region = "region", adj_info = list(), ...) {
  
  if (!validate_model(model = model)) {
    stop("Select a valid model.")
  }
  
  model <- substring(text = model, first = c(1, 4), last = c(3, 7))[2]
  
  if (model %in% c("ABST", "ABXT", "ABSS", "ABTT", "XBXT", "AASS", "AATT", "BBSS", "BBTT", "ABYZ", "ABXZ", "ABYY", "ABZZ", "XBXZ", "AAYY", "AAZZ", "BBYY", "BBZZ") & length(adj_info) == 0) {
    stop("Please, provide 'adj_info'.")
  }
  
  if (model %in% c("AACC", "AADD", "BBCC", "BBDD", "AASS", "AATT", "BBSS", "BBTT", "AAYY", "AAZZ", "BBYY", "BBZZ", "AAXX", "BBXX", "AXXX")) {
    if (!length(nonlinear) == 0) {
      stop("You cannot include non-linear effects for this model.")
    }
    m <- substring(text = model, first = c(1, 2), last = c(1, 4))[1]
    if (m == "A") {
      cov <- cov.tilde
    } 
  }
  
  # Pre-computed quantities
  observed <- which(data[[obs]] == 1)
  N <- nrow(data)
  N_obs <- length(observed)
  pre_computed <- list(observed = observed, N = N, N_obs = N_obs)
  
  if (model == "ABCD") {
    d <- data_ABCD(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region)
  } else if (model == "ABXD") {
    d <- data_ABXD(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region)
  } else if (model == "ABCC" | model == "ABDD") {
    d <- data_ABCC(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region)
  } else if (model == "XBXD") {
    d <- data_XBXD(data = data, pre_computed = pre_computed,                        cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region)
  } else if (model == "AACC" | model == "AADD" | model == "BBCC" | model == "BBDD") {
    d <- data_AACC(data = data, pre_computed = pre_computed,                        cov = cov,                                 time = time, pop.haz = pop.haz, region = region)
  } else if (model == "ABST") {
    d <- data_ABST(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "ABXT") {
    d <- data_ABXT(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "ABSS" | model == "ABTT") {
    d <- data_ABSS(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "XBXT") {
    d <- data_XBXT(data = data, pre_computed = pre_computed,                        cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "AASS" | model == "AATT" | model == "BBSS" | model == "BBTT") {
    d <- data_AASS(data = data, pre_computed = pre_computed,                        cov = cov,                                 time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "ABYZ") {
    d <- data_ABYZ(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "ABXZ") {
    d <- data_ABXZ(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "ABYY" | model == "ABZZ") {
    d <- data_ABYY(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "XBXZ") {
    d <- data_XBXZ(data = data, pre_computed = pre_computed,                        cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "AAYY" | model == "AAZZ" | model == "BBYY" | model == "BBZZ") {
    d <- data_AAYY(data = data, pre_computed = pre_computed,                        cov = cov,                                 time = time, pop.haz = pop.haz, region = region, adj_info  = adj_info)
  } else if (model == "ABXX") {
    d <- data_ABXX(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde, cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz)
  } else if (model == "XBXX") {
    d <- data_XBXX(data = data, pre_computed = pre_computed,                        cov = cov, nonlinear = nonlinear, df = df, time = time, pop.haz = pop.haz)
  } else if (model == "AAXX" | model == "BBXX") {
    d <- data_AAXX(data = data, pre_computed = pre_computed,                        cov = cov,                                 time = time, pop.haz = pop.haz)
  } else if (model == "AXXX") {
    d <- data_AXXX(data = data, pre_computed = pre_computed, cov.tilde = cov.tilde,                                            time = time, pop.haz = pop.haz)
  }
  
  d
}

data_ABCD <- function (data, pre_computed, cov.tilde, cov, nonlinear, df, time, pop.haz, region, ...) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov.tilde)
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X, X = dm$X, N_reg = length(unique(data[[region]])), region = as.integer(data[[region]]))
  
}

data_ABXD <- data_ABCD
data_ABCC <- data_ABCD

data_XBXD <- function (data, pre_computed, cov, nonlinear, df, time, pop.haz, region, ...) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X, N_reg = length(unique(data[[region]])), region = as.integer(data[[region]]))
  
}

data_AACC <- function (data, pre_computed, cov, time, pop.haz, region, ...) {
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X, N_reg = length(unique(data[[region]])), region = as.integer(data[[region]]))
  
}

data_ABST <- function (data, pre_computed, cov.tilde, cov, nonlinear, df, time, pop.haz, region, adj_info, ...) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov.tilde)
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X, X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]))
  
}

data_ABXT <- data_ABST
data_ABSS <- data_ABST

data_XBXT <- function (data, pre_computed, cov, nonlinear, df, time, pop.haz, region, adj_info, ...) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]))
  
}

data_AASS <- function (data, pre_computed, cov, time, pop.haz, region, adj_info, ...) {
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]))
  
}

data_ABYZ <- function (data, pre_computed, cov.tilde, cov, nonlinear, df, time, pop.haz, region, adj_info, ...) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov.tilde)
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X, X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]), scaling_factor = adj_info$scaling_factor)
  
}

data_ABXZ <- data_ABYZ
data_ABYY <- data_ABYZ

data_XBXZ <- function (data, pre_computed, cov, nonlinear, df, time, pop.haz, region, adj_info, ...) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]), scaling_factor = adj_info$scaling_factor)
  
}

data_AAYY <- function (data, pre_computed, cov, time, pop.haz, region, adj_info, ...) {
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X, N_reg = adj_info$N_reg, N_edges = adj_info$N_edges, node1 = adj_info$node1, node2 = adj_info$node2, region = as.integer(data[[region]]), scaling_factor = adj_info$scaling_factor)
  
}

data_ABXX <- function (data, pre_computed, cov.tilde, cov, nonlinear, df, time, pop.haz, ...) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov.tilde)
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X, X = dm$X)
  
}

data_XBXX <- function (data, pre_computed, cov, nonlinear, df, time, pop.haz) {
  
  spl <- matrix_spline(data = data, nonlinear = nonlinear, df = df)
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov, spl = spl)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, M_spl = (length(nonlinear) * df), df = df, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X)
  
}

data_AAXX <- function (data, pre_computed, cov, time, pop.haz) {
  
  dm <- design_matrix(data = data, N = pre_computed$N, cov = cov)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M = dm$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X = dm$X)
  
}

data_AXXX <- function (data, pre_computed, cov.tilde, time, pop.haz) {
  
  dm_tilde <- design_matrix(data = data, N = pre_computed$N, cov = cov.tilde)
  
  list(N = pre_computed$N, N_obs = pre_computed$N_obs, M_tilde = dm_tilde$M, obs = pre_computed$observed, time = data[[time]], pop_haz = data[[pop.haz]], X_tilde = dm_tilde$X)
  
}
