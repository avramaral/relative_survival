
# Write a function to compute S(t) based on the estimated parameters.

fit <- readRDS(file = "DATA/fitted_random_effects.rds")
fitted_data <- rstan::extract(fit)