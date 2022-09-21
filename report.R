source("header.R")

model <- "LN_ABCD"
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1])

r <- readRDS(file = paste("FITTED_MODELS/", dist, "/", model, ".rds", sep = ""))
fit <- r$fit

# Visualization

print(fit, pars = c("log_lik"), include = F)
pairs(x = fit, pars = c("log_lik", "energy__", "lp__", "v", "v_tilde", "u", "u_tilde"), include = F)
traceplot(object = fit, pars = c("log_lik", "energy__", "lp__", "v", "v_tilde", "u", "u_tilde"), include = F)

# Model comparison

models <- c("LN_ABCD", "LN_ABST", "LL_ABCD", "LL_ABST")
distributions <- sapply(X = models, FUN = function (x) { gsub(pattern = "_", replacement = "", x = substring(text = x, first = c(1, 4), last = c(3, 7))[1]) })

fits <- list()
loos <- list()

for (i in 1:length(models)) {
  print(i)
  temp <- readRDS(file = paste("FITTED_MODELS/", distributions[i], "/", models[i], ".rds", sep = ""))
  fits[[i]] <- temp$fit
  loos[[i]] <- compute_loo(fit = fits[[i]])
  # I can also include the Bayes Factor, but I have to save the "bridge_sampler()" object when fitting the model. See commit #36
}

loo_compare(loos) # Refitting the model for the outliers, or relying only on the Bayes Factor

# Result processing

data <- readRDS(file = "DATA/leuk.rds")
map <- readRDS(file = "DATA/nwengland_map.rds")

time <- seq(from = 0.025, to = 4, by = 0.025)

mu_age <- mean(data$age)
sd_age <- sd(data$age)
age <- 40 + time
age <- (age - mu_age) / sd_age

X_tilde <- matrix(data = c(age), ncol = 1, byrow = F) 
X <- matrix(data = c(age, rep(1, length(time)), rep(0.5, length(time)), rep(1.2, length(time))), ncol = 4, byrow = F) 

res <- result_processing(fit = fit, model = model, time = time, X_tilde = X_tilde, X = X)

excHaz <- res$excHaz
excCumHaz <- res$excCumHaz
netSur <- res$netSur

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
