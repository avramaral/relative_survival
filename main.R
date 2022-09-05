library(spdep)
library(rstan)
library(shinystan)
library(parallel)

source("utils.R")
source("data_stan.R")
source("distributions.R")
source("result_processing.R")

data <- readRDS(file = "DATA/data.rds")
map  <- readRDS(file = "DATA/nwengland_map.rds")

adj <- poly2nb(pl = map)
adj <- nb2mat(neighbours = adj, style = "B")
N_reg <- nrow(adj)

nodes <- adj_quantities(adj, as.numeric(rownames(adj))) # From "utils.R"
node1 <- nodes$node1
node2 <- nodes$node2

adj_info <- list(N_reg = N_reg, N_edges = length(node1), node1 = node1, node2 = node2)

model <- 6

d <- data_stan(data = data, model = model, cov_tilde = c("age"), cov = c("sex", "wbc", "dep"), intercept_tilde = T, intercept = T, adj_info = adj_info)
# d <- data_stan(data = data, intercept = T, model = model)

str(d)

### Stan Modeling

distribution <- "PGW" # PGW, LN, or LL

seed <- 1
chains <- 4
iter <- 10e3
warmup <- 8e3

start_time <- Sys.time()

fit <- stan(file = paste("MODELS/", distribution, "/", distribution, model, ".stan", sep = ""), 
            data = d,
            chains = chains,
            iter = iter,
            warmup = warmup,
            # seed = seed,
            control = list(adapt_delta = 0.80, max_treedepth = 10),
            cores = getOption(x = "mc.cores", default = detectCores())) 

end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

fitted_data <- extract(fit)

saveRDS(object = fit, file = paste("FITTED_MODELS/", distribution, "/", distribution, model, ".rds", sep = ""))
# fit <- readRDS(file = paste("FITTED_MODELS/", distribution, "/", distribution, model, ".rds", sep = ""))

# Assess Fitted Model

print(fit)
pairs(x = fit, pars = c("energy__", "lp__"), include = F) # log = T

par <- fitted_data$mu
par(family = 'LM Roman 10', mfrow = c(1, 1))
plot_chains(par = par, chains = chains, iter = iter, warmup = warmup)

launch_shinystan(fit) # Shiny App

# Result Processing

N_samples <- length(fitted_data$lp__)

time <- seq(from = 0.025, to = 4, by = 0.025)

X_tilde <- matrix(data = c(1, 1.5), nrow = length(time), ncol = 2, byrow = T) 
X <- matrix(data = c(1, 1, 0.5, 1.2), nrow = length(time), ncol = 4, byrow = T) 

res <- result_processing(model = model, fitted_data = fitted_data, N_samples = N_samples, N_reg = N_reg, distribution = distribution, time = time, X_tilde = X_tilde, X = X)

excHaz <- res$excHaz
excCumHaz <- res$excCumHaz
netSur <- res$netSur

# Result Visualization

region <- 1
par(family = 'LM Roman 10', mfrow = c(1, 1))

plot_summary_curve(time = time, obj = excHaz, region = region, ylab = "Excess Hazard", return_values = T)
plot_summary_curve(time = time, obj = excCumHaz, region = region, ylab = "Excess Cumulative Hazard", return_values = T)
plot_summary_curve(time = time, obj = netSur, region = region, return_values = T)

plot_all_regions(time = time, obj = excHaz, N_reg = N_reg, ylab = "Excess Hazard", return_values = T)
plot_all_regions(time = time, obj = excCumHaz, N_reg = N_reg, ylab = "Excess Cumulative Hazard", pos_legend = "topleft", return_values = T)
plot_all_regions(time = time, obj = netSur, N_reg = N_reg, return_values = T)
