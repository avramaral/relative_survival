
###########################
### Auxiliary Functions ###
###########################

### "main.R" Functions

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

### "data_stan.R" Functions

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

### "result_processing.R" Functions

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

add_re <- function (model, lp, lp_tilde, random_effects, ...) {
  
  u_tilde <- random_effects[1]
  u <- random_effects[2]
  
  if (model == 1) {
    lp_re_tilde <- add_re_aux(lp = lp_tilde, random_effect = u_tilde)
    lp_re <- add_re_aux(lp = lp, random_effect = u)
  } else if (model == 2) {
    lp_re_tilde <- lp_tilde
    lp_re <- add_re_aux(lp = lp, random_effect = u)
  } else if (model == 3) {
    lp_re_tilde <- add_re_aux(lp = lp_tilde, random_effect = u)
    lp_re <- add_re_aux(lp = lp, random_effect = u)
  } else if (model == 4) {
    lp_re_tilde <- lp_tilde
    lp_re <- add_re_aux(lp = lp, random_effect = u)
  } else if (model == 5) {
    lp_re_tilde <- add_re_aux(lp = lp_tilde, random_effect = u)
    lp_re <- add_re_aux(lp = lp, random_effect = u)
  } else if (model %in% c(6:9)) {
    lp_re_tilde <- lp_tilde
    lp_re <- lp
  } else {
    stop("Select a valid model.")
  }
  
  list(lp_re_tilde = lp_re_tilde, lp_re = lp_re)
}

#############################
### Plots & Visualization ###
#############################

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

plot_summary_curve <- function (time, obj, region = 1, ylab = "Net Survival", distribution = "PGW", spatial = T, return_values = F, ...) {
  
  if (!spatial) {
    region <- 1
  }
  
  M <- apply(X = obj[, , region], MARGIN = c(1), FUN = mean)
  L <- apply(X = obj[, , region], MARGIN = c(1), FUN = quantile, prob = c(0.025))
  U <- apply(X = obj[, , region], MARGIN = c(1), FUN = quantile, prob = c(0.975))
  
  plot(NA, xlim = c(0, max(time)), ylim = c(0, max(U)), xlab = "Time", ylab = ylab, main = ifelse(test = spatial, yes = paste("Region ", sprintf('%02d', region), " (", distribution, ")", sep = ""), no = paste("ALL (", distribution, ")", sep = "")))
  polygon(x = c(time, rev(time)), y = c(L, rev(U)), col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), border = F)
  lines(x = time, y = M, col = rgb(red = 0.1, green = 0.1, blue = 0.1), lty = 2)
  lines(x = time, y = L, col = rgb(red = 0.1, green = 0.1, blue = 0.1), lty = 1)
  lines(x = time, y = U, col = rgb(red = 0.1, green = 0.1, blue = 0.1), lty = 1)
  p <- recordPlot()
  
  if (return_values) { return(list(plot = p, M = M, L = L, U = U)) }
}

plot_all_regions <- function (time, obj, N_reg, ylab = "Net Survival", distribution = "PGW", spatial = T, pos_legend = "topright", return_values = F, ...) {
  
  if (!spatial) {
    N_reg <- 1
  } 
  
  meanCurves <- apply(X = obj, MARGIN = c(1, 3), FUN = mean)

  plot(NA, xlim = c(0, max(time)), ylim = c(0, max(meanCurves)), xlab = "Time", ylab = ylab, main = paste("All regions (", distribution, ")", sep = ""))
  for (j in 1:N_reg) {
    lines(time, meanCurves[, j], col = j)
  }
  if (spatial) {
    legend(x = pos_legend, inset = 0.01, legend = paste("Region ", sprintf('%02d', c(1:N_reg)), sep = ""), col = 1:N_reg, lty = 1, box.lty = 0, cex = 0.5)
  }
  p <- recordPlot()
  
  if (return_values) { return(list(plot = p, meanCurves = meanCurves)) }
}
