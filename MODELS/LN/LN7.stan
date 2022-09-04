#include ../stan_functions.stan

data {
  int<lower = 1> N;
  int<lower = 0> N_obs;
  int<lower = 0> M;
  int<lower = 1, upper = N> obs[N_obs];
  vector<lower = 0>[N] time;
  vector<lower = 0>[N] pop_haz;
  matrix[N, M] X;
}

parameters {
  vector[M] beta;
  
  real mu;
  real log_sigma;
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  {
    vector[N] lp;
    
    vector[N] excessHaz;
    vector[N] cumExcessHaz;
    
    lp = linear_predictor(N, X, beta);
    
    excessHaz = hazLN(N, time, mu, exp(log_sigma), 0) .* exp(lp);
    cumExcessHaz = cumHazLN(N, time, mu, exp(log_sigma)) .* exp(lp);
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  for (i in 1:M) { target += normal_lpdf(beta[i] | 0, 1); }
  
  // LN location parameters
  target += normal_lpdf(mu | 0, 1); 
  
  // LN scale parameters
  target += cauchy_lpdf(log_sigma | 0, 1); // Check all the priors
  
}

// generated quantities { } 
