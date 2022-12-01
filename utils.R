# "main.R" functions

adj_list <- function(map = NULL, W = NULL, sf = FALSE, ...) {
  
  if (!class(map) == "SpatialPolygons") {
    if (is.matrix(W)) {
      adj <- W  
    } else {
     stop("'W' must be a neighbor matrix.") 
    }
  } else {
    adj <- poly2nb(pl = map)
    adj <- nb2mat(neighbours = adj, style = "B") 
  }
  
  N_reg <- nrow(adj)
  rownames(adj) <- 1:N_reg
  colnames(adj) <- 1:N_reg
  
  nodes <- adj_quantities(adj, as.numeric(rownames(adj)))
  node1 <- nodes$node1
  node2 <- nodes$node2
  
  if (sf) {
    adj.sparse <- sparseMatrix(i = node1, j = node2, x = 1, symmetric = TRUE)
    Q <- Diagonal(n = N_reg, x = rowSums(adj.sparse)) - adj.sparse
    Q_pert <- Q + Diagonal(n = N_reg) * max(diag(Q)) * sqrt(.Machine$double.eps)
    Q_inv <- inla.qinv(Q = Q_pert, constr = list(A = matrix(data = 1, nrow = 1, ncol = N_reg), e = 0))
    scaling_factor <- exp(mean(log(diag(Q_inv))))
    adj_info <- list(N_reg = N_reg, N_edges = length(node1), node1 = node1, node2 = node2, scaling_factor = scaling_factor)
  } else {
    adj_info <- list(N_reg = N_reg, N_edges = length(node1), node1 = node1, node2 = node2)
  }
  adj_info
}

adj_quantities <- function (adj, unique_regions, ...) {
  node1 <- c()
  node2 <- c()
  for (i in 2:(dim(adj)[1])) { # Lower Triangular Matrix 
    for (j in 1:(i - 1)) {
      if (adj[i, j] != 0) { 
        node1 <- c(node1, unique_regions[i])
        node2 <- c(node2, unique_regions[j])
      }
    }
  }
  list(node1 = node1, node2 = node2)
}

compute_loo <- function (fit, ...) {
  log_lik <- extract_log_lik(stanfit = fit, merge_chains = F)
  r_eff <- relative_eff(exp(log_lik), cores = getOption(x = "mc.cores", default = detectCores()))
  loo <- loo(x = log_lik, r_eff = r_eff, cores = getOption(x = "mc.cores", default = detectCores()))
  loo
}

make_dummy <- function (var, name = tail(str_split(string = deparse(substitute(var)), pattern = "\\$")[[1]], 1), ...) {
  n_categories <- length(unique(var))
  m <- matrix(0, nrow = length(var), ncol = n_categories)
  
  for(i in 1:length(var)) {
    m[i, var[i]] <- 1
  }
  colnames(m) <- paste(rep(x = name, times = n_categories), "_", 1:n_categories, sep = "")
  m
}

# "data_stan.R" functions

validate_model <- function (model, ...) {
  model %in% ALL_MDL
}

matrix_spline <- function (data, nonlinear, df, ...) {
  
  X <- matrix(data = 0, nrow = nrow(data), ncol = 0)
  if (length(nonlinear)) {
    for (i in 1:length(nonlinear)) {
      X <- cbind(X, bs(x = data[[nonlinear[i]]], df = df)[, 1:df])
    }
    colnames(X) <- paste(rep(nonlinear, each = df), rep(1:df, times = length(nonlinear)), sep = "_")
  }
  
  X
}

design_matrix <- function (data, N, cov, spl = matrix(data = 0, nrow = 1, ncol = 0), ...) {
  X_names <- cov
  int <- data.frame()[1:N, ]
  X <- as.matrix(cbind(int, data[X_names]))
  rownames(X) <- NULL
  M <- ncol(X)
  if (ncol(spl)) { X <- cbind(spl, X); M <- M + ncol(spl) }
  list(X = X, M = M)
}

# "result_processing.R" Functions

check_n_cov <- function (var, ...) {
  if (is.null(var)) {
    m <- 0
  } else {
    m <- ncol(var)
  }
  m
}

compute_lp <- function (m, X, coeff, ...) {
  if (m != 0) {
    lp <- X %*% coeff
  } else {
    lp <- as.matrix(rep(x = 0, times = nrow(X)))
  }
  lp
}

add_re_aux <- function (lp, random_effect, ...) {
  if (!is.null(random_effect)) {
    res <- lp + random_effect
  } else {
    res <- lp
  }
  res
}

add_re <- function (fitted_data, model, lp, lp_tilde, i, j, ...) {
  
  if (model %in% c("ABCD", "ABXD", "ABCC", "ABDD", "XBXD", "AACC", "AADD", "BBCC", "BBDD")) {
    re_tilde <- ifelse(test = is.null(fitted_data$v_tilde[i, j]), yes = 0, no = fitted_data$v_tilde[i, j])
    re       <- ifelse(test = is.null(fitted_data$v[i, j]),       yes = 0, no = fitted_data$v[i, j])
  } else if (model %in% c("ABST", "ABXT", "ABSS", "ABTT", "XBXT", "AASS", "AATT", "BBSS", "BBTT")) {
    re_tilde <- ifelse(test = is.null(fitted_data$u_tilde[i, j]), yes = 0, no = fitted_data$u_tilde[i, j])
    re       <- ifelse(test = is.null(fitted_data$u[i, j]),       yes = 0, no = fitted_data$u[i, j])
  } else if (model %in% c("ABYZ", "ABXZ", "ABYY", "ABZZ", "XBXZ", "AAYY", "AAZZ", "BBYY", "BBZZ")) {
    re_tilde <- ifelse(test = is.null(fitted_data$convolved_re_tilde[i, j]), yes = 0, no = fitted_data$convolved_re_tilde[i, j])
    re       <- ifelse(test = is.null(fitted_data$convolved_re[i, j]),       yes = 0, no = fitted_data$convolved_re[i, j])
  } else {
    re_tilde <- 0
    re       <- 0
  }
  
  if (model %in% c("ABCC", "ABDD", "AACC", "AADD", "BBCC", "BBDD", "ABSS", "ABTT", "AASS", "AATT", "BBSS", "BBTT", "ABYY", "ABZZ", "AAYY", "AAZZ", "BBYY", "BBZZ")) {
    lp_re_tilde <- add_re_aux(lp = lp_tilde, random_effect = re)
    lp_re       <- add_re_aux(lp = lp, random_effect = re)
  } else {
    lp_re_tilde <- add_re_aux(lp = lp_tilde, random_effect = re_tilde)
    lp_re       <- add_re_aux(lp = lp, random_effect = re)
  }
  
  list(lp_re_tilde = lp_re_tilde, lp_re = lp_re)
}

add_re_mod <- function (fitted_data, model, lp, lp_tilde, i, j, ...) {
  
  re_per_ind <- function (re, region, ...) {
    res <- c()
    for (i in 1:length(region)) {
      aux <- re[region[i]]
      res <- c(res, aux)
    }
    res
  }
  
  if (model %in% c("ABCD", "ABXD", "ABCC", "ABDD", "XBXD", "AACC", "AADD", "BBCC", "BBDD")) {
    re_tilde <- re_per_ind(re = fitted_data$v_tilde[i, ], region = j)
    re       <- re_per_ind(re = fitted_data$v[i, ], region = j)
  } else if (model %in% c("ABST", "ABXT", "ABSS", "ABTT", "XBXT", "AASS", "AATT", "BBSS", "BBTT")) {
    re_tilde <- re_per_ind(re = fitted_data$u_tilde[i, ], region = j)
    re       <- re_per_ind(re = fitted_data$u[i, ], region = j)
  } else if (model %in% c("ABYZ", "ABXZ", "ABYY", "ABZZ", "XBXZ", "AAYY", "AAZZ", "BBYY", "BBZZ")) {
    re_tilde <- re_per_ind(re = fitted_data$convolved_re_tilde[i, ], region = j)
    re       <- re_per_ind(re = fitted_data$convolved_re[i, ], region = j)
  } else {
    re_tilde <- 0
    re       <- 0
  }
  
  if (model %in% c("ABCC", "ABDD", "AACC", "AADD", "BBCC", "BBDD", "ABSS", "ABTT", "AASS", "AATT", "BBSS", "BBTT", "ABYY", "ABZZ", "AAYY", "AAZZ", "BBYY", "BBZZ")) {
    lp_re_tilde <- lp_tilde + re
    lp_re       <- lp + re
  } else {
    lp_re_tilde <- lp_tilde + re_tilde
    lp_re       <- lp + re
  }
  
  list(lp_re_tilde = lp_re_tilde, lp_re = lp_re)
} 

# Visualization

plot_chains <- function (par, chains, iter, warmup) {
  n_iter <- (iter - warmup)
  plot(x = NA, xlim = c(1, n_iter), ylim = c(min(par), max(par)), xlab = "Iterations", ylab = "Parameter")
  for (i in 1:chains) {
    lines(x = par[((i - 1) * n_iter + 1):(n_iter * i)], col = i)
  }
  legend(x = "bottomright", inset = 0, legend = paste("Chain ", sprintf('%02d', c(1:chains)), sep = ""), col = 1:chains, lty = 1, box.lty = 1, cex = 0.75)
  p <- recordPlot()
  p
}

plot_summary_curve <- function (time, obj, region = 1, ylab = "Net Survival", dist = "PGW", return_values = F, ...) {
  
  M <- apply(X = obj[, , region], MARGIN = c(1), FUN = mean)
  L <- apply(X = obj[, , region], MARGIN = c(1), FUN = quantile, prob = c(0.025))
  U <- apply(X = obj[, , region], MARGIN = c(1), FUN = quantile, prob = c(0.975))
  
  plot(NA, xlim = c(0, max(time)), ylim = c(0, max(U)), xlab = "Time", ylab = ylab, main = ifelse(test = dim(obj)[3] != 1, yes = paste("Region ", sprintf('%02d', region), " (", dist, ")", sep = ""), no = paste("ALL (", dist, ")", sep = "")))
  polygon(x = c(time, rev(time)), y = c(L, rev(U)), col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), border = F)
  lines(x = time, y = M, col = rgb(red = 0.1, green = 0.1, blue = 0.1), lty = 2)
  lines(x = time, y = L, col = rgb(red = 0.1, green = 0.1, blue = 0.1), lty = 1)
  lines(x = time, y = U, col = rgb(red = 0.1, green = 0.1, blue = 0.1), lty = 1)
  p <- recordPlot()
  
  if (return_values) { return(list(plot = p, M = M, L = L, U = U)) }
}

plot_all_regions <- function (time, obj, N_reg, ylab = "Net Survival", dist = "PGW", pos_legend = "topright", return_values = F, ...) {
  
  N_reg <- dim(obj)[3]
  
  meanCurves <- apply(X = obj, MARGIN = c(1, 3), FUN = mean)
  
  plot(NA, xlim = c(0, max(time)), ylim = c(0, max(meanCurves)), xlab = "Time", ylab = ylab, main = paste("All regions (", dist, ")", sep = ""))
  for (j in 1:N_reg) {
    lines(time, meanCurves[, j], col = j)
  }
  if (N_reg > 1) {
    legend(x = pos_legend, inset = 0.01, legend = paste("Region ", sprintf('%02d', c(1:N_reg)), sep = ""), col = 1:N_reg, lty = 1, box.lty = 0, cex = 0.5)
  }
  p <- recordPlot()
  
  if (return_values) { return(list(plot = p, meanCurves = meanCurves)) }
}

plot_map <- function (map, obj, t, summary = "M", title = "Net Survival", commom_legend = F, ...) {
  
  N_reg <- length(map)
  
  M <- c()
  L <- c()
  U <- c()
  
  N_fit_reg <- dim(obj)[3]
  
  for (i in 1:N_reg) {
    M <- c(M, ifelse(test = N_fit_reg == 1, yes = mean(obj[t, , 1]), no = mean(obj[t, , i])))
    L <- c(L, ifelse(test = N_fit_reg == 1, yes = quantile(obj[t, , 1], prob = c(0.025)), no = quantile(obj[t, , i], prob = c(0.025))))
    U <- c(U, ifelse(test = N_fit_reg == 1, yes = quantile(obj[t, , 1], prob = c(0.975)), no = quantile(obj[t, , i], prob = c(0.975))))
  }
  
  if (summary == "M") {
    value <- M
  } else if (summary == "L") {
    value <- L
  } else if (summary == "U") {
    value <- U
  } else {
    stop("Choose a proper summary measure.")
  }
  
  map_df <- st_as_sf(map)
  map_df$id <- 1:N_reg
  map_df$value <- value
  
  if (commom_legend) {
    min_legend <- min(c(M, L, U)) * 0.975
    max_legend <- max(c(M, L, U)) * 1.025
  } else {
    min_legend <- min(map_df$value) * 0.975
    max_legend <- max(map_df$value) * 1.025
  }
  
  p <- ggplot(data = map_df) + 
    geom_sf(aes(fill = value), color = "gray", size = 0.25) + 
    scale_fill_viridis(limits = c(min_legend, max_legend), name = title, option = "inferno") + 
    theme_bw() + 
    labs(x = "Longitude", y = "Latitude", title = title) + 
    theme(text = element_text(family = 'LM Roman 10'))
  
  p
}
