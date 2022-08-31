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

linear_predictor <- function (X, beta, ...) {
  X %*% beta
}

linear_predictor_re <- function (N, X, beta, region, u, ...) {
  res <- c()
  for (i in 1:N) {
    res <- c(res, u[as.integer(region[i])])
  }
  X %*% beta + res
}