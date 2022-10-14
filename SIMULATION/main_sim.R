source("header.R")
library(zoo)

model_data <- "LN_ABST"
sample_size <- 500
prop <- 50
load(paste("SIMULATION/DATA/", model_data, "_n_", sample_size, "_prop_", prop, ".RData", sep = ""))

N_sim <- 10

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
  r[[k]] <- fit_stan(mod = m, data = d, chains = 4, iter = 2e3, warmup = 1e3, max_treedepth = 12)
  print(r[[k]]$fit, pars = c("log_lik"), include = F)
}

saveRDS(object = r, file = paste("SIMULATION/FITTED_MODELS/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))
r <- readRDS(file = paste("SIMULATION/FITTED_MODELS/data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))

##################################################
# Plotted curves with real and estimated paramet # 
##################################################

N_reg <- adj_info$N_reg

time <- seq(from = 0.025, to = 5, by = 0.025)

mu_age <- mean(c(X_mal$age, X_fem$age))
sd_age <- sd(c(X_mal$age, X_fem$age))
age <- 72 + time
age <- (age - mu_age) / sd_age

X_tilde <- matrix(data = c(age), ncol = 1, byrow = F) 
X <- matrix(data = c(age, rep(1, length(time)), rep(0, length(time))), ncol = 3, byrow = F) # Male individuals who are 72+ yo and have mean deprivation levels

# Estimated parameters

res <- list()

excHaz    <- list()
excCumHaz <- list()
netSur    <- list()

for (k in 1:N_sim) {
  print(paste(sprintf('%03d', k), " out of ", N_sim, sep = ""))
  
  res[[k]] <- result_processing(fit = r[[k]]$fit, model = model, time = time, X_tilde = X_tilde, X = X)
  
  excHaz[[k]]    <- res[[k]]$excHaz
  excCumHaz[[k]] <- res[[k]]$excCumHaz
  netSur[[k]]    <- res[[k]]$netSur
}

compute_summary <- function (obj, ...) {
  
  M <- apply(X = obj, MARGIN = c(1), FUN = mean)
  L <- apply(X = obj, MARGIN = c(1), FUN = quantile, prob = c(0.025))
  U <- apply(X = obj, MARGIN = c(1), FUN = quantile, prob = c(0.975))
  
  list(M = M, L = L, U = U)
}

excHazALL    <- excHazSummary    <- list()
excCumHazALL <- excCumHazSummary <- list()
netSurALL    <- netSurSummary    <- list()

for (i in 1:N_reg) {
  print(paste(sprintf('%03d', i), " out of ", N_reg, sep = ""))
  for (k in 1:N_sim) {
    if (k == 1) {
      excHazALL[[i]]    <- excHaz[[k]][, , i] 
      excCumHazALL[[i]] <- excCumHaz[[k]][, , i]
      netSurALL[[i]]    <- netSur[[k]][, , i] 
    } else {
      excHazALL[[i]]    <- cbind(excHazALL[[i]], excHaz[[k]][, , i])
      excCumHazALL[[i]] <- cbind(excCumHazALL[[i]], excCumHaz[[k]][, , i])
      netSurALL[[i]]    <- cbind(netSurALL[[i]], netSur[[k]][, , i])
    }
  }
  
  excHazSummary[[i]]    <- compute_summary(excHazALL[[i]])
  excCumHazSummary[[i]] <- compute_summary(excCumHazALL[[i]])
  netSurSummary[[i]]    <- compute_summary(netSurALL[[i]])
}

# Real parameters

true_res <- true_result_processing(alpha = alpha, beta = beta, re_tilde = re_tilde, re = re, model = model_data, pars = pars, time = time, X_tilde = X_tilde, X = X)

true_excHaz <- true_res$excHaz
true_excCumHaz <- true_res$excCumHaz
true_netSur <- true_res$netSur

true_excHaz_val    <- plot_all_regions(time = time, obj = true_excHaz, N_reg = N_reg, ylab = "Excess Hazard", dist = dist, return_values = T)
true_excCumHaz_val <- plot_all_regions(time = time, obj = true_excCumHaz, N_reg = N_reg, ylab = "Excess Cumulative Hazard", dist = dist, pos_legend = "topleft", return_values = T)
true_netSur_val    <- plot_all_regions(time = time, obj = true_netSur, N_reg = N_reg, dist = dist, return_values = T)

# Manual plot comparing real and estimated curves 

for (i in 1:N_reg) {
  # esti <- excHazSummary[[i]]
  # true <- true_excHaz_val$meanCurves
  # ylab <- "Excess Hazard"
  
  esti <- netSurSummary[[i]]
  true <- true_netSur_val$meanCurves
  ylab <- "Net Survival"

  plot(NA, xlim = c(0, max(time)), ylim = c(0, max(esti$U)), xlab = "Time", ylab = ylab, main = paste("Region ", sprintf('%03d', i), sep = ""))
  polygon(x = c(time, rev(time)), y = c(esti$L, rev(esti$U)), col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), border = F)
  lines(x = time, y = esti$M, col = rgb(red = 0, green = 0, blue = 0), lty = 1)
  lines(x = time, y = true[, i], col = rgb(red = 1, green = 0, blue = 0), lty = 1)
  
  legend(x = "topright", inset = 0.01, legend = c("Estimated Curves", "True Curves"), col = c(1, 2), lty = c(1, 1), box.lty = 0)
}

# Compute difference between areas for the estimated and real curves
# As in https://stackoverflow.com/questions/4954507/calculate-the-area-under-a-curve

diff_areas <- c()
for (i in 1:N_reg) {
  # real <- true_excHaz_val$meanCurves[, i]
  # esti <- excHazSummary[[i]]$M[id]  
  
  real <- true_netSur_val$meanCurves[, i]
  esti <- netSurSummary[[i]]$M[id]
    
  id <- order(time)
  area_real <- sum(diff(time[id]) * rollmean(real, 2))
  area_esti <- sum(diff(time[id]) * rollmean(esti, 2))
  diff_areas <- c(diff_areas, abs(area_real - area_esti))
}

summary_obj <- list(excHazSummary = excHazSummary, excCumHazSummary = excCumHazSummary, netSurSummary = netSurSummary, true_excHaz_val = true_excHaz_val, true_excCumHaz_val = true_excCumHaz_val, true_netSur_val = true_netSur_val, diff_areas = diff_areas)
saveRDS(object = summary_obj, file = paste("SIMULATION/FITTED_MODELS/SUMMARY_data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))
summary_obj <- readRDS(file = paste("SIMULATION/FITTED_MODELS/SUMMARY_data_", model_data, "_fit_", model, "_n_", sample_size, "_prop_", prop, "_fit_result.rds", sep = ""))

# plot(NA, xlim = c(0, max(time)), ylim = c(0, ceiling(max(true, esti))), xlab = "Time", ylab = "Values", main = paste("Net Survival (", N_reg, " regions)", sep = ""))
# for (i in 1:N_reg) {
#   lines(time, true[, i], col = i)
#   if ((ncol(esti) != 1) | (i == 1))
#   lines(time, esti[, i], col = i, lty = 2)
# }
# legend(x = "topright", inset = 0.01, legend = c("True Curves", "Estimated Curves"), col = c(1, 1), lty = c(1, 2), box.lty = 0)



