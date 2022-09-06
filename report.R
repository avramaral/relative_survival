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

plot_all_curves <- function (time, res, N_reg, distribution, ...) {
  excHaz    <- plot_all_regions(time = time, obj = res$excHaz,    N_reg = N_reg, ylab = "Excess Hazard", distribution = distribution, return_values = T)
  excCumHaz <- plot_all_regions(time = time, obj = res$excCumHaz, N_reg = N_reg, ylab = "Excess Cumulative Hazard", distribution = distribution, pos_legend = "topleft", return_values = T)
  netSur    <- plot_all_regions(time = time, obj = res$netSur,    N_reg = N_reg, distribution = distribution, return_values = T)
  
  list(excHaz = excHaz, excCumHaz = excCumHaz, netSur = netSur)
}

### START HERE

model <- 6
N_reg <- 24
spatial <- ifelse(test = N_reg <= 5, yes = T, no = F)

###

fit_PGW <- readRDS(file = paste("FITTED_MODELS/PGW/PGW", model, ".rds", sep = "")); fitted_data_PGW <- extract(fit_PGW)
fit_LN  <- readRDS(file = paste("FITTED_MODELS/LN/LN"  , model, ".rds", sep = "")); fitted_data_LN  <- extract(fit_LN )
fit_LL  <- readRDS(file = paste("FITTED_MODELS/LL/LL"  , model, ".rds", sep = "")); fitted_data_LL  <- extract(fit_LL )

print(fit_PGW, pars = c("log_lik"), include = F)
print(fit_LN , pars = c("log_lik"), include = F)
print(fit_LL , pars = c("log_lik"), include = F)

# `loo`
(loo_PGW <- compute_loo(stanfit = fit_PGW))
(loo_LN  <- compute_loo(stanfit = fit_LN ))
(loo_LL  <- compute_loo(stanfit = fit_LL ))
(loo_result <- loo_compare(list(PGW = loo_PGW, LN = loo_LN, LL = loo_LL)))

# Bridge for 'Bayes factor'
bridge_PGW <- readRDS(file = paste("FITTED_MODELS/PGW/bridge_PGW", model, ".rds", sep = ""))
bridge_LN  <- readRDS(file = paste("FITTED_MODELS/LN/bridge_LN"  , model, ".rds", sep = ""))
bridge_LL  <- readRDS(file = paste("FITTED_MODELS/LL/bridge_LL"  , model, ".rds", sep = ""))

(bf_PGW_LN <- bayes_factor(x1 = bridge_PGW, x2 = bridge_LN))
(bf_PGW_LL <- bayes_factor(x1 = bridge_PGW, x2 = bridge_LL))
(bf_LN_LL  <- bayes_factor(x1 = bridge_LN , x2 = bridge_LL))

### VISUALIZATION

par(family = 'LM Roman 10', mfrow = c(1, 1))

N_samples <- length(fitted_data_PGW$lp__)

time <- seq(from = 0.025, to = 4, by = 0.025)

X_tilde <- matrix(data = c(1.5), nrow = length(time), ncol = 1, byrow = T) 
X <- matrix(data = c(1.5, 1, 0.5, 1.2), nrow = length(time), ncol = 4, byrow = T) 

res_PGW <- result_processing(model = model, fitted_data = fitted_data_PGW, N_samples = N_samples, N_reg = N_reg, distribution = "PGW", time = time, X_tilde = X_tilde, X = X)
res_LN  <- result_processing(model = model, fitted_data = fitted_data_LN , N_samples = N_samples, N_reg = N_reg, distribution = "LN" , time = time, X_tilde = X_tilde, X = X)
res_LL  <- result_processing(model = model, fitted_data = fitted_data_LL , N_samples = N_samples, N_reg = N_reg, distribution = "LL" , time = time, X_tilde = X_tilde, X = X)

res_plot_PGW <- plot_ind_curves(res = res_PGW, N_reg = N_reg, time = time, distribution = "PGW", spatial = spatial)
res_plot_LN  <- plot_ind_curves(res = res_LN , N_reg = N_reg, time = time, distribution = "LN" , spatial = spatial)
res_plot_LL  <- plot_ind_curves(res = res_LL , N_reg = N_reg, time = time, distribution = "LL" , spatial = spatial)

if (spatial) {
  all_curves_PGW <- plot_all_curves(time = time, res = res_PGW, N_reg = N_reg, distribution = "PGW")
  all_curves_LN  <- plot_all_curves(time = time, res = res_LN , N_reg = N_reg, distribution = "LN" )
  all_curves_LL  <- plot_all_curves(time = time, res = res_LL , N_reg = N_reg, distribution = "LL" )
}

# Draw maps
