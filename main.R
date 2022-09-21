source("header.R")

data <- readRDS(file = "DATA/leuk.rds")

# Optional
data$age <- scale(data$age)
data$wbc <- scale(data$wbc)
data$dep <- scale(data$dep)

map <- readRDS(file = "DATA/nwengland_map.rds")
adj_info <- adj_list(map = map)

model <- "PGWAACC"
dist <- gsub(pattern = "_", replacement = "", x = substring(text = model, first = c(1, 4), last = c(3, 7))[1])

d <- data_stan(data = data, model = model, cov.tilde = c("age"), cov = c("age", "wbc", "sex", "dep"), nonlinear = c(), adj_info = adj_info)
r <- fit_stan(data = d, model = model)

print(r$fit, pars = c("log_lik"), include = F)
saveRDS(object = r, file = paste("FITTED_MODELS/", dist, "/", model, ".rds", sep = ""))

r$model
r$time_taken



# Bayes Factor
bridge <- bridge_sampler(samples = r$fit, cores = getOption(x = "mc.cores", default = detectCores()), silent = T)
saveRDS(object = bridge, file = paste("FITTED_MODELS/", dist, "/bridge_", model, ".rds", sep = ""))
