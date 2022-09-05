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
  
  vector[N_reg] u_tilde;
  vector[N_reg] u;
}

transformed parameters {
    vector[N] lp_tilde;
    vector[N] lp;
  
    vector[N] excessHaz;
    vector[N] cumExcessHaz;
    
    lp_tilde = linear_predictor_re(N, X_tilde, alpha, region, u_tilde);
    lp = linear_predictor_re(N, X, beta, region, u);
    
    excessHaz = hazLN(N, time .* exp(lp_tilde), mu, exp(log_sigma), 0) .* exp(lp);
    cumExcessHaz = cumHazLN(N, time .* exp(lp_tilde), mu, exp(log_sigma)) .* exp(lp - lp_tilde);
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
  
  // LN locations parameters
  target += normal_lpdf(mu | 0, 1); 
  
  // LN scale parameters
  target += cauchy_lpdf(log_sigma | 0, 1); // Check all the priors
  
  // Random effects
  target += icar_normal_lpdf(u_tilde | N_reg, node1, node2);
  target += icar_normal_lpdf(u | N_reg, node1, node2);
  
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
