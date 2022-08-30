library(spdep)
library(rstan)
library(parallel)

source("utils.R")
source("result_processing.R")

data <- readRDS(file = "DATA/data.rds")
map  <- readRDS(file = "DATA/nwengland_map.rds")

adj <- poly2nb(pl = map)
adj <- nb2mat(neighbours = adj, style = "B")
N_reg <- nrow(adj)

nodes <- adj_quantities(adj, as.numeric(rownames(adj))) # From "utils.R"
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
                  N_reg = N_reg,
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

# saveRDS(object = fit, file = "DATA/fitted_random_effects.rds")
# fit <- readRDS(file = "DATA/fitted_random_effects.rds")

fitted_data <- rstan::extract(fit)

### For Generated Quantities

N_gen <- 101
range_time <- range(data$time)
new_time <- seq(from = floor(range_time[1]), to = ceiling(range_time[2]), length.out = N_gen)

new_age <- 1 # after scaling it
new_sex <- 0
new_wbc <- 2
new_dep <- 1
X_tilde_gen <- as.matrix(cbind(1, new_age))
X_gen <- as.matrix(cbind(1, new_sex, new_wbc, new_dep))

N_samples <- (iter - warmup) * chains

lp_tilde <- array(data = 0, dim = c(N_samples, N_reg))
lp <- array(data = 0, dim = c(N_samples, N_reg))
estHaz <- array(data = 0, dim = c(N_gen, N_samples, N_reg))
estCumHaz <- array(data = 0, dim = c(N_gen, N_samples, N_reg))
estSurviv <- array(data = 0, dim = c(N_gen, N_samples, N_reg))

progressbar <- txtProgressBar(min = 1, max = N_reg, initial = 1)
for (i in 1:N_reg) {
  for (j in 1:N_samples) {
    lp_tilde[j, i] <- linear_predictor(X = X_tilde_gen, fixed_coeff = fitted_data$alpha[j, ], region = as.integer(i), random_effect = fitted_data$u_tilde[j, ])
    lp[j, i] <- linear_predictor(X = X_gen, fixed_coeff = fitted_data$beta[j, ], region = as.integer(i), random_effect = fitted_data$u[j, ])
    estHaz[, j, i] <- hazPGW(time = (new_time * exp(lp_tilde[j, i])), eta = fitted_data$eta[j], nu = fitted_data$nu[j], theta = fitted_data$theta[j]) * exp(lp[j, i])
    estCumHaz[, j, i] <- cumHazPGW(time = (new_time * exp(lp_tilde[j, i])), eta = fitted_data$eta[j], nu = fitted_data$nu[j], theta = fitted_data$theta[j]) * exp(lp[j, i] - lp_tilde[j, i])
    estSurviv[, j, i] <- estSurv(cumHaz = estCumHaz[, j, i])
  }
  setTxtProgressBar(progressbar, i)
}
close(progressbar)

# Excess Survival plot for a specific region (still missing to include the population hazard)
new_region <- 10
par(family = 'LM Roman 10', mfrow = c(1, 1))
plot(x = NA, y = NA, main = paste("Region ",  sprintf('%02d', new_region), sep = ""), xlab = "Time", ylab = "Excess Survival", xlim = c(new_time[1], tail(new_time, 1)), ylim = c(0, 1))
for (j in 6001:8000) { # Samples were arbitrarily chosen as the fourth chain and region "new_region"
  lines(x = new_time, y = estSurviv[, j, new_region], col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.2))
}

# Excess Survival plot for all regions with averaged curves
estSurviv_curves <- apply(X = estSurviv[, 6001:8000, ], MARGIN = c(1, 3), FUN = mean) # Samples were arbitrarily chosen as the fourth chain
par(family = 'LM Roman 10', mfrow = c(1, 1))
plot(x = NA, y = NA, main = paste("All regions"), xlab = "Time", ylab = "Excess Survival", xlim = c(new_time[1], tail(new_time, 1)), ylim = c(0, 1))
for (i in (1:ncol(estSurviv_curves))) {
  lines(x = new_time, y = estSurviv_curves[, i], col = i)
}
