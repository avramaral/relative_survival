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

model {
  // --------------
  // Log-likelihood
  // --------------
  
  {
    vector[N] lp;
    
    vector[N] excessHaz;
    vector[N] cumExcessHaz;
    
    lp = linear_predictor(N, X, beta);
    
    excessHaz = hazPGW(N, time .* exp(lp), exp(log_eta), exp(log_nu), exp(log_theta), 0) .* exp(lp);
    cumExcessHaz = cumHazPGW(N, time .* exp(lp), exp(log_eta), exp(log_nu), exp(log_theta));
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  for (i in 1:M) { target += normal_lpdf(beta[i] | 0, 1); }
  
  // PGW scale parameters
  target += cauchy_lpdf(log_eta | 0, 1); 
  
  // PGW shape parameters
  target += cauchy_lpdf(log_nu | 0, 1);
  target += cauchy_lpdf(log_theta | 0, 1); // Check all the priors
  
}

// generated quantities { } 
