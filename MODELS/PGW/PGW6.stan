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

parameters {
  vector[M_tilde] alpha;
  vector[M] beta;
  
  real log_eta;
  real log_nu;
  real<lower = 0> theta;
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  {
    vector[N] lp_tilde;
    vector[N] lp;
  
    vector[N] excessHaz;
    vector[N] cumExcessHaz;
    
    lp_tilde = linear_predictor(N, X_tilde, alpha);
    lp = linear_predictor(N, X, beta);
    
    excessHaz = hazPGW(N, time .* exp(lp_tilde), exp(log_eta), exp(log_nu), theta, 0) .* exp(lp);
    cumExcessHaz = cumHazPGW(N, time .* exp(lp_tilde), exp(log_eta), exp(log_nu), theta) .* exp(lp - lp_tilde);
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  for (i in 1:M_tilde) { target += normal_lpdf(alpha[i] | 0, 1); }
  for (i in 1:M) { target += normal_lpdf(beta[i] | 0, 1); }
  
  // PGW scale parameters
  target += cauchy_lpdf(log_eta | 0, 1); 
  
  // PGW shape parameters
  target += cauchy_lpdf(log_nu | 0, 1);
  target += gamma_lpdf(theta | 0.75, 0.75); // Check all the priors
  
}

// generated quantities { }
