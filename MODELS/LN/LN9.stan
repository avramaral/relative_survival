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
  
  real mu;
  real log_sigma;
}

transformed parameters {
  
  real<lower=0> sigma;
  
  sigma = exp(log_sigma);

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
    
    excessHaz = hazLN(N, time .* exp(lp_tilde), mu, sigma, 0);
    cumExcessHaz = cumHazLN(N, time .* exp(lp_tilde), mu, sigma) .* exp(-lp_tilde);
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  alpha ~ normal(0, 10);
  
  // LN location parameters
  target += normal_lpdf(mu | 0, 10); 
  
  // LN scale parameters
  target += cauchy_lpdf(log_sigma | 0, 1); // Check all the priors
  
}

// generated quantities { } 
