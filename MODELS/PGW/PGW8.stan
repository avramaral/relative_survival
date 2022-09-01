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
  
  real log_eta;
  real log_nu;
  real log_theta;
}

transformed parameters {
  
  real<lower=0> eta;
  real<lower=0> nu; 
  real<lower=0> theta;
  
  eta = exp(log_eta);
  nu = exp(log_nu);
  theta = exp(log_theta);
  
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
    
    excessHaz = hazPGW(N, time .* exp(lp), eta, nu, theta, 0) .* exp(lp);
    cumExcessHaz = cumHazPGW(N, time .* exp(lp), eta, nu, theta);
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  beta ~ normal(0, 100);
  
  // PGW scale parameters
  target += cauchy_lpdf(log_eta | 0, 5); 
  
  // PGW shape parameters
  target += cauchy_lpdf(log_nu | 0, 5);
  target += cauchy_lpdf(log_theta | 0, 5); // Check all the priors
  
}

// generated quantities { } 
