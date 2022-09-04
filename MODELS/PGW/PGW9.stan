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

parameters {
  vector[M_tilde] alpha;
  
  real log_eta;
  real log_nu;
  real log_theta;
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  {
    vector[N] lp_tilde;
    
    vector[N] excessHaz;
    vector[N] cumExcessHaz;
    
    lp_tilde = linear_predictor(N, X_tilde, alpha);
    
    excessHaz = hazPGW(N, time .* exp(lp_tilde), exp(log_eta), exp(log_nu), exp(log_theta), 0);
    cumExcessHaz = cumHazPGW(N, time .* exp(lp_tilde), exp(log_eta), exp(log_nu), exp(log_theta)) .* exp(-lp_tilde);
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  for (i in 1:M_tilde) { target += normal_lpdf(alpha[i] | 0, 1); }
  
  // PGW scale parameters
  target += cauchy_lpdf(log_eta | 0, 1); 
  
  // PGW shape parameters
  target += cauchy_lpdf(log_nu | 0, 1);
  target += cauchy_lpdf(log_theta | 0, 1); // Check all the priors
  
}

// generated quantities { } 
