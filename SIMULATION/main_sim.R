args <- commandArgs(trailingOnly = TRUE)
model_data  <- args[1]
model       <- args[2]
prop        <- args[3]
sample_size <- args[4]
N_sim_fit   <- args[5]

library("methods")
source("header.R")
library("zoo")

manual_effects <- FALSE    #####

# model_data  <- "LN_ABST" #####
# sample_size <- 5000      #####
# prop        <- 75        #####

if (!manual_effects) {
  load(paste("SIMULATION/DATA/", model_data, "_n_", sample_size, "_prop_", prop, ".RData", sep = ""))
} else {
  load(paste("SIMULATION/DATA/SPECIAL/", model_data, "_n_", sample_size, "_prop_", prop, ".RData", sep = ""))
}

N_sim <- as.numeric(N_sim_fit) #####
N_valid_fit <- 100             #####
min_N_valid_fit_N_sim <- min(N_valid_fit, N_sim)

##################################################

adj_list_simplified <- function(W, ...) {
  adj <- W
  N_reg <- nrow(adj)
  
  nodes <- adj_quantities(adj, 1:9)
  node1 <- nodes$node1
  node2 <- nodes$node2
  
  list(N_reg = N_reg, N_edges = length(node1), node1 = node1, node2 = node2)
}

W <- rbind(c(0, 1, 1, 0, 0, 0, 0, 0, 0),
           c(1, 0, 1, 1, 1, 0, 0, 0, 0),
           c(1, 1, 0, 1, 0, 0, 0, 0, 0),
           c(0, 1, 1, 0, 1, 1, 0, 1, 0),
           c(0, 1, 0, 1, 0, 0, 0, 1, 1),
           c(0, 0, 0, 1, 0, 0, 1, 1, 0),
           c(0, 0, 0, 0, 0, 1, 0, 1, 0),
           c(0, 0, 0, 1, 1, 1, 1, 0, 1),
           c(0, 0, 0, 0, 1, 0, 0, 1, 0))

adj_info <- adj_list_simplified(W)

# model <- "LN_ABCD" #####
dist    <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1])

if (dist %in% c("LN", "LL")) {
  analyzed_pars <- c("mu", "sigma", "alpha", "beta")
} else if (dist %in% c("PGW", "GG")) {
  analyzed_pars <- c("eta", "nu", "theta", "alpha", "beta")
} else if (dist == "GAM") {
  analyzed_pars <- c("eta", "nu", "alpha", "beta")
}

m <- compile_model(model = model)

r <- list()
bridge <- list()
count <- 0
for (k in 1:N_sim) {
  set.seed((999 + k))
  
  print(paste("Current valid fitted models: ", count, sep = ""))
  if (count < min_N_valid_fit_N_sim) {
    print(paste(sprintf('%03d', k), " out of ", N_sim, sep = ""))
    d <- data_stan(data = data[[k]], model = model, cov.tilde = c("age"), cov = c("age", "sex", "dep"), nonlinear = c(), adj_info = adj_info)
    r_temp <- fit_stan(mod = m, data = d, chains = 4, iter = 4e3, warmup = 2e3, max_treedepth = 12, adapt_delta = 0.8)
    
    if (as.logical(sum(summary(r_temp$fit, pars = analyzed_pars)$summary[, "Rhat"] > 1.5))) {
      valid_fit <- FALSE
    } else {
      valid_fit <- TRUE
    }
    
    if (valid_fit) {
      r[[(count + 1)]] <- r_temp
      bridge[[(count + 1)]] <- bridge_sampler(samples = r[[(count + 1)]]$fit, cores = getOption(x = "mc.cores", default = detectCores()), silent = T)
      count <- count + 1
    }
  }
}

if (!manual_effects) {
  saveRDS(object = r, file = paste("SIMULATION/FITTED_MODELS/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))
  saveRDS(object = bridge, file = paste("SIMULATION/FITTED_MODELS/BRIDGE/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_bridge.rds", sep = ""))
} else {
  saveRDS(object = r, file = paste("SIMULATION/FITTED_MODELS/SPECIAL/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))
  saveRDS(object = bridge, file = paste("SIMULATION/FITTED_MODELS/SPECIAL/BRIDGE/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_bridge.rds", sep = ""))
}

print(paste("Number of valid fitted models: ", count, sep = ""))
