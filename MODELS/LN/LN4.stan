#include ../stan_functions.stan

data {
  int<lower = 1> N;
  int<lower = 0> N_obs;
  int<lower = 0> M;
  int<lower = 1, upper = N> obs[N_obs];
  vector<lower = 0>[N] time;
  vector<lower = 0>[N] pop_haz;
  matrix[N, M] X;
  
  // Information about the adjacency matrix
  int<lower = 0> N_reg;
  int<lower = 0> N_edges;
  int<lower = 1, upper = N_reg> node1[N_edges]; // Extremes of edges
  int<lower = 1, upper = N_reg> node2[N_edges];
  int<lower = 1, upper = N_reg> region[N]; // Region for each observation
}

parameters {
  vector[M] beta;
  
  real mu;
  real log_sigma;

  vector[N_reg] u;
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
    
    lp = linear_predictor_re(N, X, beta, region, u);
    
    excessHaz = hazLN(N, time, mu, sigma, 0) .* exp(lp);
    cumExcessHaz = cumHazLN(N, time, mu, sigma) .* exp(lp);
    
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
  
  // Random effects
  target += icar_normal_lpdf(u | N_reg, node1, node2);
  
}

// generated quantities { } 
