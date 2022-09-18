#include ../stan_functions.stan

data {
  int<lower = 1> N;
  int<lower = 0> N_obs;
  int<lower = 0> M_tilde;
  int<lower = 0> M;
  int<lower = 0> M_spl;
  int<lower = 0> df;
  int<lower = 1, upper = N> obs[N_obs];
  vector<lower = 0>[N] time;
  vector<lower = 0>[N] pop_haz;
  matrix[N, M_tilde] X_tilde;
  matrix[N, M] X;
  
  int<lower = 0> N_reg;
  int<lower = 1, upper = N_reg> region[N]; 
}

transformed data {
  vector[N] ind_obs;
  ind_obs = rep_vector(0.0, N);
  ind_obs[obs] = rep_vector(1.0, N_obs);
  
  int<lower = 0> N_spl;
  N_spl = M_spl %/% df;
  
  array[N_spl] matrix[N, M_spl] X_spl; 
  array[N_spl] matrix[M_spl, M_spl] B; 
  
  if (N_spl != 0) {
    for (i in 1:N_spl) {
      X_spl[i] = X[, ((i - 1) * df + 1):(i * df)];
      B[i] = ((N - (0.5 *  (N - N_obs))) / df) * inverse_spd(X_spl[i]' * X_spl[i]); 
    }
  }
}

parameters {
  vector[M_tilde] alpha;
  vector[M] beta;
  
  real<lower = 0> eta;
  real<lower = 0> nu;

  vector[N_reg] v;
  
  real log_sigma_v;
  
  vector[N_spl] log_sigma_B;
}

transformed parameters {
  vector[N] lp_tilde;
  vector[N] lp;
  
  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  lp_tilde = linear_predictor(N, X_tilde, alpha);
  lp = linear_predictor_re(N, X, beta, region, v);
  
  excessHaz = hazGAM(N, time .* exp(lp_tilde), eta, nu, 0) .* exp(lp);
  cumExcessHaz = cumHazGAM(N, time .* exp(lp_tilde), eta, nu) .* exp(lp - lp_tilde);
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  target += sum(log(pop_haz[obs] + excessHaz[obs])) - sum(cumExcessHaz);
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Non-linear fixed coefficients
  if (N_spl != 0) {
    for (i in 1:N_spl) {
      target += multi_normal_lpdf(beta[((i - 1) * df + 1):(i * df)] | rep_vector(0.0, df), pow(exp(log_sigma_B[i]), 2) * B[i]);
    }
  }

  // Linear Fixed coefficients
  for (i in 1:M_tilde) { 
    target += normal_lpdf(alpha[i] | 0, 10); 
  }
  
  for (i in (M_spl + 1):M) { 
    target += normal_lpdf(beta[i] | 0, 10); 
  }
  
  // GAM scale parameters
  target += cauchy_lpdf(eta | 0, 1); 
  
  // GAM shape parameters
  target += cauchy_lpdf(nu | 0, 1);
  
  // Random effects
  target += normal_lpdf(v | 0, 1);
  
  // Hyperpriors
  target += normal_lpdf(log_sigma_v | 0, 1); 
  
  if (N_spl != 0) {
    for (i in 1:N_spl) {
      target += normal_lpdf(log_sigma_B[i] | 0, 1);
    }
  }
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
