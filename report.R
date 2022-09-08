source("header.R")

compute_loo <- function (stanfit, ...) {
  log_lik <- extract_log_lik(stanfit = stanfit, merge_chains = F)
  r_eff <- relative_eff(exp(log_lik), cores = getOption(x = "mc.cores", default = detectCores()))
  loo <- loo(x = log_lik, r_eff = r_eff, cores = getOption(x = "mc.cores", default = detectCores()))
  loo
}

plot_ind_curves <- function (res, N_reg, time, distribution, spatial = T, ...) {
  
  if (!spatial) {
    N_reg <- 1
  } 
  
  res_excHaz <- list(); res_excCumHaz <- list(); res_netSur <- list();
  
  for (i in 1:N_reg) {
    res_excHaz[[i]]    <- plot_summary_curve(time = time, obj = res$excHaz,    region = i, ylab = "Excess Hazard", distribution = distribution, spatial = spatial, return_values = T)
    res_excCumHaz[[i]] <- plot_summary_curve(time = time, obj = res$excCumHaz, region = i, ylab = "Excess Cumulative Hazard", distribution = distribution, spatial = spatial, return_values = T)
    res_netSur[[i]]    <- plot_summary_curve(time = time, obj = res$netSur,    region = i, distribution = distribution, spatial = spatial, return_values = T)
  }
  
  list(res_excHaz = res_excHaz, res_excCumHaz = res_excCumHaz, res_netSur = res_netSur)
}

plot_all_curves <- function (time, res, N_reg, distribution, spatial = T, ...) {
  excHaz    <- plot_all_regions(time = time, obj = res$excHaz,    N_reg = N_reg, ylab = "Excess Hazard", distribution = distribution, spatial = spatial, return_values = T)
  excCumHaz <- plot_all_regions(time = time, obj = res$excCumHaz, N_reg = N_reg, ylab = "Excess Cumulative Hazard", distribution = distribution, spatial = spatial, pos_legend = "topleft", return_values = T)
  netSur    <- plot_all_regions(time = time, obj = res$netSur,    N_reg = N_reg, distribution = distribution, spatial = spatial, return_values = T)
  
  list(excHaz = excHaz, excCumHaz = excCumHaz, netSur = netSur)
}

### START HERE

mu_age <- 60.72579
sd_age <- 18.33427

model <- 1
N_reg <- 24
spatial <- ifelse(test = model <= 5, yes = T, no = F)

fit_PGW <- list(); fit_LN <- list(); fit_LL <- list();
fitted_data_PGW <- list(); fitted_data_LN <- list(); fitted_data_LL <- list();
loo_PGW <- list(); loo_LN <- list(); loo_LL <- list();
bridge_PGW <- list(); bridge_LN <- list(); bridge_LL <- list()
res_PGW <- list(); res_LN <- list(); res_LL <- list()
res_plot_PGW <- list(); res_plot_LN <- list(); res_plot_LL <- list()

###

fit_PGW[[model]] <- readRDS(file = paste("FITTED_MODELS/PGW/PGW", model, ".rds", sep = "")); fitted_data_PGW[[model]] <- extract(fit_PGW[[model]])
fit_LN[[model]]  <- readRDS(file = paste("FITTED_MODELS/LN/LN"  , model, ".rds", sep = "")); fitted_data_LN[[model]]  <- extract(fit_LN[[model]] )
fit_LL[[model]]  <- readRDS(file = paste("FITTED_MODELS/LL/LL"  , model, ".rds", sep = "")); fitted_data_LL[[model]]  <- extract(fit_LL[[model]] )

print(fit_PGW[[model]], pars = c("log_lik"), include = F)
print(fit_LN[[model]] , pars = c("log_lik"), include = F)
print(fit_LL[[model]] , pars = c("log_lik"), include = F)

# `loo`
(loo_PGW[[model]] <- compute_loo(stanfit = fit_PGW[[model]]))
(loo_LN[[model]]  <- compute_loo(stanfit = fit_LN[[model]] ))
(loo_LL[[model]]  <- compute_loo(stanfit = fit_LL[[model]] ))
loo_compare(list(PGW = loo_PGW[[model]], LN = loo_LN[[model]], LL = loo_LL[[model]])) # Include the desired models for comparison
# loo_compare(list(LL_1 = loo_LL[[1]], LL_6 = loo_LL[[6]]))

# Bridge for 'Bayes factor'
bridge_PGW[[model]] <- readRDS(file = paste("FITTED_MODELS/PGW/bridge_PGW", model, ".rds", sep = ""))
bridge_LN[[model]]  <- readRDS(file = paste("FITTED_MODELS/LN/bridge_LN"  , model, ".rds", sep = ""))
bridge_LL[[model]]  <- readRDS(file = paste("FITTED_MODELS/LL/bridge_LL"  , model, ".rds", sep = ""))

bayes_factor(x1 = bridge_LL[[model]], x2 = bridge_LN[[model]] )
bayes_factor(x1 = bridge_LL[[model]], x2 = bridge_PGW[[model]])
bayes_factor(x1 = bridge_LN[[model]], x2 = bridge_PGW[[model]]) # Include the desired models for comparison
# bayes_factor(x1 = bridge_LL[[1]], x2 = bridge_LL[[6]])

### VISUALIZATION

par(family = 'LM Roman 10', mfrow = c(1, 1))

N_samples <- length(fitted_data_PGW[[model]]$lp__)

time <- seq(from = 0.025, to = 4, by = 0.025)

age <- 40 + time
age <- (age - mu_age) / sd_age

X_tilde <- matrix(data = c(age), ncol = 1, byrow = F) 
X <- matrix(data = c(age, rep(1, length(time)), rep(0.5, length(time)), rep(1.2, length(time))), ncol = 4, byrow = F) 

res_PGW[[model]] <- result_processing(model = model, fitted_data = fitted_data_PGW[[model]], N_samples = N_samples, N_reg = N_reg, distribution = "PGW", time = time, X_tilde = X_tilde, X = X, spatial = spatial)
res_LN[[model]]  <- result_processing(model = model, fitted_data = fitted_data_LN[[model]] , N_samples = N_samples, N_reg = N_reg, distribution = "LN" , time = time, X_tilde = X_tilde, X = X, spatial = spatial)
res_LL[[model]]  <- result_processing(model = model, fitted_data = fitted_data_LL[[model]] , N_samples = N_samples, N_reg = N_reg, distribution = "LL" , time = time, X_tilde = X_tilde, X = X, spatial = spatial)

res_plot_PGW[[model]] <- plot_ind_curves(res = res_PGW[[model]], N_reg = N_reg, time = time, distribution = "PGW", spatial = spatial)
res_plot_LN[[model]]  <- plot_ind_curves(res = res_LN[[model]] , N_reg = N_reg, time = time, distribution = "LN" , spatial = spatial)
res_plot_LL[[model]]  <- plot_ind_curves(res = res_LL[[model]] , N_reg = N_reg, time = time, distribution = "LL" , spatial = spatial)
# print(res_plot_PGW[[1]]$res_netSur[[24]]$plot)

if (spatial) {
  plot_all_curves(time = time, res = res_PGW[[model]], N_reg = N_reg, distribution = "PGW", spatial = spatial)
  plot_all_curves(time = time, res = res_LN[[model]] , N_reg = N_reg, distribution = "LN" , spatial = spatial)
  plot_all_curves(time = time, res = res_LL[[model]] , N_reg = N_reg, distribution = "LL" , spatial = spatial)
}

# Maps

library(sp)
library(ggplot2)
library(viridis)

plot_map <- function (map, values, title = "Net Survival", ...) {
  
  N_reg <- length(map)
  
  map_df <- st_as_sf(map)
  map_df$id <- 1:N_reg
  map_df$value <- values
  
  p <- ggplot(data = map_df) + 
       geom_sf(aes(fill = value), color = "gray", size = 0.25) + 
       scale_fill_viridis(limits = c(min(map_df$value) * 0.95, max(map_df$value) * 1.05), name = title, option = "inferno") + 
       theme_bw() + 
       labs(x = "Longitude", y = "Latitude", title = title) + 
       theme(text = element_text(family = 'LM Roman 10'))
  
  p
}

map  <- readRDS(file = "DATA/nwengland_map.rds")

t <- which(time == 1)
values <- c()
for (i in 1:N_reg) {
  values <- c(values, res_plot_PGW[[model]]$res_netSur[[i]]$M[t])
  # values <- colMeans(fitted_data_PGW[[model]]$u) # Random Effects (or "u_tilde")
}

plot_map(map = map, values = values, title = "Net Survival")
