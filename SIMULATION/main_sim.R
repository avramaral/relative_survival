source("header.R")
library(zoo)

model_data <- "LN_ABCD"
sample_size <- 2000
prop <- 50
load(paste("SIMULATION/DATA/", model_data, "_n_", sample_size, "_prop_", prop, ".RData", sep = ""))

N_sim <- 1

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

model <- "LN_ABCD"
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1])

m <- compile_model(model = model)

r <- list()
for (k in 1:N_sim) {
  print(paste(sprintf('%03d', k), " out of ", N_sim, sep = ""))
  d <- data_stan(data = data[[k]], model = model, cov.tilde = c("age"), cov = c("age", "sex", "dep"), nonlinear = c(), adj_info = adj_info)
  r[[k]] <- fit_stan(mod = m, data = d, chains = 4, iter = 2e3, warmup = 1e3, max_treedepth = 10)
  print(r[[k]]$fit, pars = c("log_lik"), include = F)
}

saveRDS(object = r, file = paste("SIMULATION/FITTED_MODELS/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))
# r <- readRDS(file = paste("SIMULATION/FITTED_MODELS/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))
