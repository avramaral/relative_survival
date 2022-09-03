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

transformed parameters {
  
  real<lower=0> sigma;
  
  sigma = exp(log_sigma);

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
    
    excessHaz = hazLN(N, time .* exp(lp), mu, sigma, 0) .* exp(lp);
    cumExcessHaz = cumHazLN(N, time .* exp(lp), mu, sigma);
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  beta ~ normal(0, 10);
  
  // LN location parameters
  target += normal_lpdf(mu | 0, 10); 
  
  // LN scale parameters
  target += cauchy_lpdf(log_sigma | 0, 1); // Check all the priors
  
}

// generated quantities { } 
