#include ../stan_functions.stan

data {
  int<lower = 1> N;
  int<lower = 0> N_obs;
  int<lower = 0> M_tilde;
  int<lower = 1, upper = N> obs[N_obs];
  vector<lower = 0>[N] time;
  vector<lower = 0>[N] pop_haz;
  matrix[N, M_tilde] X_tilde;
}

transformed data {
  vector[N] ind_obs;
  ind_obs = rep_vector(0.0, N);
  ind_obs[obs] = rep_vector(1.0, N_obs);
}

parameters {
  vector[M_tilde] alpha;
  
  real mu;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[N] lp_tilde;
  
  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  lp_tilde = linear_predictor(N, X_tilde, alpha);
  
  excessHaz = hazLN(N, time .* exp(lp_tilde), mu, sigma, 0);
  cumExcessHaz = cumHazLN(N, time .* exp(lp_tilde), mu, sigma) .* exp(-lp_tilde);
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Linear Fixed coefficients
  for (i in 1:M_tilde) { 
    target += normal_lpdf(alpha[i] | 0, 10); 
  }
  
  // LN location parameters
  target += normal_lpdf(mu | 0, 10); 
  
  // LN scale parameters
  target += cauchy_lpdf(sigma | 0, 1);
  
}

generated quantities { 
  vector[N] log_lik;
  for (i in 1:N) {
    if (ind_obs[i] == 0.0) {
      log_lik[i] = - cumExcessHaz[i];
    } else {
      log_lik[i] = log(pop_haz[i] + excessHaz[i]) - cumExcessHaz[i];
    }
  }
}
