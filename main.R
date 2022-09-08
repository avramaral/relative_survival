source("header.R")

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
spatial <- ifelse(test = model <= 5, yes = T, no = F)

d <- data_stan(data = data, model = model, cov_tilde = c("age"), cov = c("age", "sex", "wbc", "dep"), adj_info = adj_info)
# d <- data_stan(data = data, model = model)

str(d)

### Stan Modeling

distribution <- "PGW" # PGW, LN, or LL

seed <- 1
chains <- 4
iter <- 4e3
warmup <- 2e3

start_time <- Sys.time()

fit <- stan(file = paste("MODELS/", distribution, "/", distribution, model, ".stan", sep = ""), 
            data = d,
            chains = chains,
            iter = iter,
            warmup = warmup,
            # seed = seed,
            pars = c("lp_tilde", "lp", "excessHaz", "cumExcessHaz"),
            include = F,
            control = list(adapt_delta = 0.80, max_treedepth = 10),
            cores = getOption(x = "mc.cores", default = detectCores())) 

end_time <- Sys.time()
time_taken <- end_time - start_time
time_taken

fitted_data <- extract(fit)

saveRDS(object = fit, file = paste("FITTED_MODELS/", distribution, "/", distribution, model, ".rds", sep = ""))
# fit <- readRDS(file = paste("FITTED_MODELS/", distribution, "/", distribution, model, ".rds", sep = ""))

### Assess Fitted Model

print(fit, pars = c("log_lik"), include = F)
pairs(x = fit, pars = c("log_lik", "energy__", "lp__"), include = F) # log = T

par <- fitted_data$mu
par(family = 'LM Roman 10', mfrow = c(1, 1))
plot_chains(par = par, chains = chains, iter = iter, warmup = warmup)

## Model Comparison

# "loo"
log_lik <- extract_log_lik(stanfit = fit, merge_chains = F)
r_eff <- relative_eff(exp(log_lik), cores = getOption(x = "mc.cores", default = detectCores()))
loo <- loo(x = log_lik, r_eff = r_eff, cores = getOption(x = "mc.cores", default = detectCores()))
print(loo)
# loo_compare(list("M1" = loo1, "M2" = loo2))

# "Bayes factor"
bridge <- bridge_sampler(samples = fit, cores = getOption(x = "mc.cores", default = detectCores()), silent = T)
saveRDS(object = bridge, file = paste("FITTED_MODELS/", distribution, "/bridge_", distribution, model, ".rds", sep = ""))
# bf <- bayes_factor(x1 = bridge1, x2 = bridge2)

### Result Processing

N_samples <- length(fitted_data$lp__)

time <- seq(from = 0.025, to = 4, by = 0.025)

mu_age <- attr(x = data$age, which = "scaled:center")
sd_age <- attr(x = data$age, which = "scaled:scale")
age <- 40 + time
age <- (age - mu_age) / sd_age

X_tilde <- matrix(data = c(age), ncol = 1, byrow = F) 
X <- matrix(data = c(age, rep(1, length(time)), rep(0.5, length(time)), rep(1.2, length(time))), ncol = 4, byrow = F) 

res <- result_processing(model = model, fitted_data = fitted_data, N_samples = N_samples, N_reg = N_reg, distribution = distribution, time = time, X_tilde = X_tilde, X = X, spatial = spatial)

excHaz <- res$excHaz
excCumHaz <- res$excCumHaz
netSur <- res$netSur

### Result Visualization

region <- 1
par(family = 'LM Roman 10', mfrow = c(1, 1))

plot_summary_curve(time = time, obj = excHaz, region = region, ylab = "Excess Hazard", distribution = distribution, spatial = spatial, return_values = T)
plot_summary_curve(time = time, obj = excCumHaz, region = region, ylab = "Excess Cumulative Hazard", distribution = distribution, spatial = spatial, return_values = T)
plot_summary_curve(time = time, obj = netSur, region = region, distribution = distribution, spatial = spatial, return_values = T)

if (spatial) {
  plot_all_regions(time = time, obj = excHaz, N_reg = N_reg, ylab = "Excess Hazard", distribution = distribution, spatial = spatial, return_values = T)
  plot_all_regions(time = time, obj = excCumHaz, N_reg = N_reg, ylab = "Excess Cumulative Hazard", distribution = distribution, spatial = spatial, pos_legend = "topleft", return_values = T)
  plot_all_regions(time = time, obj = netSur, N_reg = N_reg, distribution = distribution, spatial = spatial, return_values = T)
}
