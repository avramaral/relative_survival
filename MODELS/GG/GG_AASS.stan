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
  int<lower = 0> N_edges;
  int<lower = 1, upper = N_reg> node1[N_edges]; 
  int<lower = 1, upper = N_reg> node2[N_edges];
  int<lower = 1, upper = N_reg> region[N]; 
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

  vector[N_reg] u;
}

transformed parameters {
  vector[N] lp;
  
  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  lp = linear_predictor_re(N, X, beta, region, u);
  
  excessHaz = hazGG(N, time .* exp(lp), eta, nu, theta, 0) .* exp(lp);
  cumExcessHaz = cumHazGG(N, time .* exp(lp), eta, nu, theta);
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
  
  // GG scale parameters
  target += cauchy_lpdf(eta | 0, 1); 
  
  // GG shape parameters
  target += cauchy_lpdf(nu | 0, 1);
  target += gamma_lpdf(theta | 0.65, 1 / 1.83); 
  
  // Random effects
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
