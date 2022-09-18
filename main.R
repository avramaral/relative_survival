source("header.R")

data <- readRDS(file = "DATA/leuk.rds")

# Optional
data$age <- scale(data$age)
data$wbc <- scale(data$wbc)
data$dep <- scale(data$dep)

map <- readRDS(file = "DATA/nwengland_map.rds")
adj_info <- adj_list(map = map)

model <- "PGWABCD"
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1])

d <- data_stan(data = data, model = model, cov.tilde = c("age"), cov = c("age", "wbc", "sex", "dep"), nonlinear = c(), adj_info = adj_info)
r <- fit_stan(data = d, model = model)

fit <- r$fit
print(fit, pars = c("log_lik"), include = F)

# saveRDS(object = fit, file = paste("FITTED_MODELS/", dist, "/", model, ".rds", sep = ""))

## Model comparison
loo <- compute_loo(fit = fit)
# loo_compare(list(M1 = loo1, M2 = loo2))

bridge <- bridge_sampler(samples = fit, cores = getOption(x = "mc.cores", default = detectCores()), silent = T)
# bayes_factor(x1 = bridge1, x2 = bridge2)
# saveRDS(object = bridge, file = paste("FITTED_MODELS/", dist, "/bridge_", model, ".rds", sep = ""))

## Result processing

time <- seq(from = 0.025, to = 4, by = 0.025)

mu_age <- attr(x = data$age, which = "scaled:center")
sd_age <- attr(x = data$age, which = "scaled:scale")
age <- 40 + time
age <- (age - mu_age) / sd_age

X_tilde <- matrix(data = c(age), ncol = 1, byrow = F) 
X <- matrix(data = c(age, rep(1, length(time)), rep(0.5, length(time)), rep(1.2, length(time))), ncol = 4, byrow = F) 

res <- result_processing(fit = fit, model = model, time = time, X_tilde = X_tilde, X = X)

excHaz    <- res$excHaz
excCumHaz <- res$excCumHaz
netSur    <- res$netSur

## Visualization

region <- 1
par(family = 'LM Roman 10', mfrow = c(1, 1))

plot_summary_curve(time = time, obj = excHaz, region = region, ylab = "Excess Hazard", dist = dist, return_values = T)
plot_summary_curve(time = time, obj = excCumHaz, region = region, ylab = "Excess Cumulative Hazard", dist = dist, return_values = T)
plot_summary_curve(time = time, obj = netSur, region = region, dist = dist, return_values = T)

plot_all_regions(time = time, obj = excHaz, N_reg = N_reg, ylab = "Excess Hazard", dist = dist, return_values = T)
plot_all_regions(time = time, obj = excCumHaz, N_reg = N_reg, ylab = "Excess Cumulative Hazard", dist = dist, pos_legend = "topleft", return_values = T)
plot_all_regions(time = time, obj = netSur, N_reg = N_reg, dist = dist, return_values = T)

par(family = 'LM Roman 10', mfrow = c(1, 3))
pL <- plot_map(map = map, obj = netSur, t = which(time == 1), summary = "L", title = "Net Survival (.025)", commom_legend = T)
pM <- plot_map(map = map, obj = netSur, t = which(time == 1), summary = "M", title = "Net Survival (mean)", commom_legend = T)
pU <- plot_map(map = map, obj = netSur, t = which(time == 1), summary = "U", title = "Net Survival (.975)", commom_legend = T)
grid.arrange(pL, pM, pU, nrow = 1, ncol = 3)
