#include ../stan_functions.stan

data {
  int<lower = 1> N;
  int<lower = 0> N_obs;
  int<lower = 0> M;
  int<lower = 1, upper = N> obs[N_obs];
  vector<lower = 0>[N] time;
  vector<lower = 0>[N] pop_haz;
  matrix[N, M] X;
  
  int<lower = 0> N_reg;
  int<lower = 1, upper = N_reg> region[N]; 
}

transformed data {
  vector[N] ind_obs;
  ind_obs = rep_vector(0.0, N);
  ind_obs[obs] = rep_vector(1.0, N_obs);
}

parameters {
  vector[M] beta;
  
  real log_eta;
  real log_nu;
  real<lower = 0> theta;

  vector[N_reg] v;
  
  real log_sigma_v;
}

transformed parameters {
  vector[N] lp;
  
  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  lp = linear_predictor_re(N, X, beta, region, v);
  
  excessHaz = hazPGW(N, time .* exp(lp), exp(log_eta), exp(log_nu), theta, 0) .* exp(lp);
  cumExcessHaz = cumHazPGW(N, time .* exp(lp), exp(log_eta), exp(log_nu), theta);
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
  target += cauchy_lpdf(log_eta | 0, 1); 
  
  // PGW shape parameters
  target += cauchy_lpdf(log_nu | 0, 1);
  target += gamma_lpdf(theta | 0.65, 1 / 1.83); 
  
  // Random effects
  target += normal_lpdf(v | 0, exp(log_sigma_v));
  
  // Hyperpriors
  target += normal_lpdf(log_sigma_v | 0, 1); 
  
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
