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

transformed data {
  vector[N] ind_obs;
  ind_obs = rep_vector(0.0, N);
  ind_obs[obs] = rep_vector(1.0, N_obs);
}

parameters {
  vector[M] beta;
  
  real<lower = 0> eta;
  real<lower = 0> nu;
  real<lower = 0> theta;
}

transformed parameters {
  vector[N] lp;
  
  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  lp = linear_predictor(N, X, beta);
  
  excessHaz = hazPGW(N, time .* exp(lp), eta, nu, theta, 0) .* exp(lp);
  cumExcessHaz = cumHazPGW(N, time .* exp(lp), eta, nu, theta);
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Linear Fixed coefficients
  for (i in 1:M) { 
    target += normal_lpdf(beta[i] | 0, 10); 
  }
  
  // PGW scale parameters
  target += cauchy_lpdf(eta | 0, 1); 
  
  // PGW shape parameters
  target += cauchy_lpdf(nu | 0, 1);
  target += gamma_lpdf(theta | 0.65, 1 / 1.83); 
  
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
