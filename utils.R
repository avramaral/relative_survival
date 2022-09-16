
# "main.R" functions

adj_list <- function(map, ...) {
  
  if (!class(map) == "SpatialPolygons") {
    stop("'map' is not a 'SpatialPolygons' object.")
  }
  
  adj <- poly2nb(pl = map)
  adj <- nb2mat(neighbours = adj, style = "B")
  N_reg <- nrow(adj)
  
  nodes <- adj_quantities(adj, as.numeric(rownames(adj)))
  node1 <- nodes$node1
  node2 <- nodes$node2
  
  adj_info <- list(N_reg = N_reg, N_edges = length(node1), node1 = node1, node2 = node2)
  
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