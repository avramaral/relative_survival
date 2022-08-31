
# IGNORE THIS FOR NOW

fitted_data <- rstan::extract(fit)

fitted_data





time <- seq(from = 0, to = 10, by = 0.1)
N <- length(time) 

age <- 1.5
sex <- 0
wbc <- 0.25
dep <- 0.5

M <- 5
X <- matrix(data = c(1, age, sex, wbc, dep), nrow = N, ncol = M, byrow = T)
beta <- fitted_data$

linear_predictor_re(N = N, X = X, beta = )
  
do.call(paste("cumHaz", distribution, sep = ""), N = N, time = )






get(paste("cumHaz", distribution, sep = ""))
