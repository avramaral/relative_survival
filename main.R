library(spdep)
library(rstan)
library(parallel)

source("utils.R")
source("data_stan.R")

data <- readRDS(file = "DATA/data.rds")
map  <- readRDS(file = "DATA/nwengland_map.rds")

adj <- poly2nb(pl = map)
adj <- nb2mat(neighbours = adj, style = "B")
N_reg <- nrow(adj)

nodes <- adj_quantities(adj, as.numeric(rownames(adj))) # From "utils.R"
node1 <- nodes$node1
node2 <- nodes$node2

adj_info <- list(N_reg = N_reg, N_edges = length(node1), node1 = node1, node2 = node2)

model <- 1

d <- data_stan(data = data, model = model, cov_tilde = c("age"), cov = c("sex", "wbc", "dep"), intercept_tilde = T, intercept = T, adj_info = adj_info)
str(d)

### Stan Modeling

distribution <- "PGW"

seed <- 1
chains <- 4
iter <- 20e3
warmup <- 18e3

start_time <- Sys.time()

fit <- stan(file = paste("MODELS/", distribution, "/", distribution, model, ".stan", sep = ""), 
            data = d,
            chains = chains,
            iter = iter,
            warmup = warmup,
            # seed = seed,
            control = list(adapt_delta = 0.80),
            cores = getOption(x = "mc.cores", default = detectCores())) 

end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

# saveRDS(object = fit, file = paste("FITTED_MODELS/", distribution, "/", distribution, model, ".rds", sep = ""))
# fit <- readRDS(file = paste("FITTED_MODELS/", distribution, "/", distribution, model, ".rds", sep = ""))

fitted_data <- extract(fit)
N_samples <- length(fitted_data$lp__)

check_n_cov <- function (var) {
  if (is.null(var)) {
    m <- 0
  } else {
    m <- ncol(var)
  }
  m
}

lp <- function (m, int, cov, coeff) {
  if (m != 0) {
    if (int) {
      cov <- c(1, cov)
    }
    lp <- cov %*% coeff
  } else {
    lp <- matrix(0, 1, 0)
  }
  as.numeric(lp)
}

# re <- function (model) {
#   if (model == 1) {
#     
#   } else if (model == 2)
# } 

# Set it as empty initially in the function
cov_tilde <- c(1.5) # "age"
cov <- c(1, 0.5, 1.2) # "sex", "wbc", "dep"
intercept <- T
intercept_tilde <- T

m_tilde <- check_n_cov(var = fitted_data$alpha)
m <- check_n_cov(var = fitted_data$beta)

if (((length(cov_tilde) + as.integer(intercept_tilde)) != m_tilde) | ((length(cov) + as.integer(intercept)) != m)) {
  stop("Provide the correct values for covariates.")
}

for (i in 1:N_samples) {
  lp_tilde <- compute_design_mat(m = m_tilde, int = intercept_tilde, cov = cov_tilde, coeff = fitted_data$alpha[i, ])
  lp <- compute_design_mat(m = m, int = intercept, cov = cov, coeff = fitted_data$beta[i, ])
  
  # NOW I HAVE TO INCLUDE THE RANDOM EFFECTS SOMEHOW
}





















