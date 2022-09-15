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
  real log_sigma;
}

transformed parameters {
  vector[N] lp_tilde;
  
  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  lp_tilde = linear_predictor(N, X_tilde, alpha);
  
  excessHaz = hazLL(N, time .* exp(lp_tilde), mu, exp(log_sigma), 0);
  cumExcessHaz = cumHazLL(N, time .* exp(lp_tilde), mu, exp(log_sigma)) .* exp(-lp_tilde);
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  for (i in 1:M_tilde) { target += normal_lpdf(alpha[i] | 0, 10); }
  
  // LL location parameters
  target += normal_lpdf(mu | 0, 10); 
  
  // LL scale parameters
  target += normal_lpdf(log_sigma | 0, 1); // Check all the priors
  
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
