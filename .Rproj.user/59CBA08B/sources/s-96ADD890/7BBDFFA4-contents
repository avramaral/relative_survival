library(spdep)
library(rstan)
library(parallel)

data <- readRDS(file = "DATA/data.rds")
map  <- readRDS(file = "DATA/nwengland_map.rds")

adj <- poly2nb(pl = map)
adj <- nb2mat(neighbours = adj, style = "B")

adj_quantities <- function (adj, unique_regions) {
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

nodes <- adj_quantities(adj, as.numeric(rownames(adj)))
node1 <- nodes$node1
node2 <- nodes$node2

### Data Preparation

cens <- which(data$cens == 1)

N <- nrow(data)
N_cens <- length(cens)

# Design matrices
# X tilde 
X_tilde_names <- c("age") # Select proper covariates
X_tilde <- as.matrix(cbind(rep(1, N), data[X_tilde_names]))
colnames(X_tilde) <- c("int", X_tilde_names)
M_tilde <- ncol(X_tilde)
# X
X_names <- c("sex", "wbc", "dep") # Excluding "age"
X <- as.matrix(cbind(rep(1, N), data[X_names]))
colnames(X) <- c("int", X_names) # Intercept + covariates
M <- ncol(X)

# For Generated Quantities
N_gen <- 101
range_time <- range(data$time)
new_t <- seq(from = floor(range_time[1]), to = ceiling(range_time[2]), length.out = N_gen)

# Stan data object
data_stan <- list(N = N,
                  N_cens = N_cens,
                  M_tilde = M_tilde,
                  M = M,
                  cens = cens,
                  time = data$time,
                  pop_haz = data$pop.haz,
                  X_tilde = X_tilde,
                  X = X,
                  N_reg = nrow(adj),
                  N_edges = length(node1),
                  node1 = node1,
                  node2 = node2,
                  region = as.integer(data$region))

### Stan Modeling

seed <- 1
chains <- 4
iter <- 20e3 
warmup <- 18e3

start_time <- Sys.time()

fit <- stan(file = "model.stan", 
            data = data_stan,
            chains = chains,
            iter = iter,
            warmup = warmup,
            seed = seed,
            control = list(adapt_delta = 0.99),
            cores = getOption(x = "mc.cores", default = detectCores())) 

end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

saveRDS(object = fit, file = "DATA/fitted_random_effects.rds")

fitted_data <- rstan::extract(fit)
