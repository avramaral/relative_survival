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

d <- data_stan(data = data, model = 1, adj_info = adj_info)
str(d)

### Stan Modeling

model_name <- "PGW1"

seed <- 1
chains <- 4
iter <- 10e3 
warmup <- 8e3

start_time <- Sys.time()

fit <- stan(file = paste(model_name, ".stan", sep = ""), 
            data = d,
            chains = chains,
            iter = iter,
            warmup = warmup,
            # seed = seed,
            control = list(adapt_delta = 0.99),
            cores = getOption(x = "mc.cores", default = detectCores())) 

end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

saveRDS(object = fit, file = paste("DATA/", model_name ,".rds", sep = ""))
# fit <- readRDS(file = paste("DATA/", model_name ,".rds", sep = ""))

fitted_data <- rstan::extract(fit)


