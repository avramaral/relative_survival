source("header.R")

data <- readRDS(file = "DATA/leuk.rds")

# Optional
data$age <- scale(data$age)
data$wbc <- scale(data$wbc)
data$dep <- scale(data$dep)

map <- readRDS(file = "DATA/nwengland_map.rds")
adj_info <- adj_list(map = map)

model <- "PGWABXX"
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1])

d <- data_stan(data = data, model = model, cov.tilde = c("age"), cov = c("wbc", "sex", "dep"), nonlinear = c(), adj_info = adj_info)
r <- fit_stan(data = d, model = model)

fit <- r$fit
print(fit, pars = c("log_lik"), include = F)

bridge <- bridge_sampler(samples = fit, cores = getOption(x = "mc.cores", default = detectCores()), silent = T)

# saveRDS(object = bridge, file = paste("FITTED_MODELS/", dist, "/bridge_", model, ".rds", sep = ""))
# saveRDS(object = fit, file = paste("FITTED_MODELS/", dist, "/", model, ".rds", sep = ""))





log_lik <- extract_log_lik(stanfit = fit, merge_chains = F)
r_eff <- relative_eff(exp(log_lik), cores = getOption(x = "mc.cores", default = detectCores()))
loo <- loo(x = log_lik, r_eff = r_eff, cores = getOption(x = "mc.cores", default = detectCores()))
print(loo)
