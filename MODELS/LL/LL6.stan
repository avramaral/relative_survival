#include ../stan_functions.stan

data {
  int<lower = 1> N;
  int<lower = 0> N_obs;
  int<lower = 0> M_tilde;
  int<lower = 0> M;
  int<lower = 1, upper = N> obs[N_obs];
  vector<lower = 0>[N] time;
  vector<lower = 0>[N] pop_haz;
  matrix[N, M_tilde] X_tilde;
  matrix[N, M] X;
}

transformed data {
  vector[N] ind_obs;
  ind_obs = rep_vector(0.0, N);
  ind_obs[obs] = rep_vector(1.0, N_obs);
}

parameters {
  vector[M_tilde] alpha;
  vector[M] beta;
  
  real mu;
  real log_sigma;
}

transformed parameters {
  vector[N] lp_tilde;
  vector[N] lp;

  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  lp_tilde = linear_predictor(N, X_tilde, alpha);
  lp = linear_predictor(N, X, beta);
  
  excessHaz = hazLL(N, time .* exp(lp_tilde), mu, exp(log_sigma), 0) .* exp(lp);
  cumExcessHaz = cumHazLL(N, time .* exp(lp_tilde), mu, exp(log_sigma)) .* exp(lp - lp_tilde);
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
  for (i in 1:M_tilde) { target += normal_lpdf(alpha[i] | 0, 1); }
  for (i in 1:M) { target += normal_lpdf(beta[i] | 0, 1); }
  
  // LL location parameters
  target += normal_lpdf(mu | 0, 1); 
  
  // LL scale parameters
  target += cauchy_lpdf(log_sigma | 0, 1); // Check all the priors
  
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
