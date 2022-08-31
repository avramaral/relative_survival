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
  
  // Information about the adjacency matrix
  int<lower = 0> N_reg;
  int<lower = 0> N_edges;
  int<lower = 1, upper = N_reg> node1[N_edges]; // Extremes of edges
  int<lower = 1, upper = N_reg> node2[N_edges];
  int<lower = 1, upper = N_reg> region[N]; // Region for each observation
}

parameters {
  vector[M_tilde] alpha;
  vector[M] beta;
  
  real log_eta;
  real log_nu;
  real log_theta;

  vector[N_reg] u;
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
    vector[N] lp_tilde;
    vector[N] lp;
    
    vector[N] excessHaz;
    vector[N] cumExcessHaz;
    
    lp_tilde = linear_predictor(N, X_tilde, alpha);
    lp = linear_predictor_re(N, X, beta, region, u);
    
    excessHaz = hazPGW(N, time .* exp(lp_tilde), eta, nu, theta, 0) .* exp(lp);
    cumExcessHaz = cumHazPGW(N, time .* exp(lp_tilde), eta, nu, theta) .* exp(lp - lp_tilde);
    
    target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  // PGW scale parameters
  target += cauchy_lpdf(log_eta | 0, 2.5); 
  
  // PGW shape parameters
  target += cauchy_lpdf(log_nu | 0, 2.5);
  target += cauchy_lpdf(log_theta | 0, 2.5); // Check all the priors
  
  // Random effects
  target += icar_normal_lpdf(u | N_reg, node1, node2);
  
}

// generated quantities { } 
