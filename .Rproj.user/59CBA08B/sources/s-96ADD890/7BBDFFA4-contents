library(spdep)
library(rstan)

data <- readRDS(file = "DATA/data.rds")
map  <- readRDS(file = "DATA/nwengland_map.rds")

adj <- poly2nb(pl = map)
adj <- nb2mat(neighbours = adj, style = "B")

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

# Stan data object
data_stan <- list(N = N,
                  N_cens = N_cens,
                  M_tilde = M_tilde,
                  M = M,
                  cens = cens,
                  time = data$time,
                  pop_haz = data$pop.haz,
                  X_tilde = X_tilde,
                  X = X #, still missing information about the adjacency matrix
                 )

### Stan Modeling

chains <- 2
iter <- 2e3
warmup <- 1e3

start.time <- Sys.time()

fit <- stan(file = "model.stan", 
            data = data_stan,
            chains = chains,
            iter = iter,
            warmup = warmup,
            # control = list(adapt_delta = 0.99)
            ) # Include parallelization

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

saveRDS(object = fit, file = "DATA/fitted_no_random_effects.rds")

fitted_data <- rstan::extract(fit)
