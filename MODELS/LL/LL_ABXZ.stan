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
  int<lower = 0> N_edges;
  int<lower = 1, upper = N_reg> node1[N_edges];
  int<lower = 1, upper = N_reg> node2[N_edges];
  int<lower = 1, upper = N_reg> region[N];
  
  real<lower = 0> scaling_factor;
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
  
  real mu;
  real<lower = 0> sigma;
  
  real<lower = 0, upper = 1> rho;

  vector[N_reg] v;
  
  vector[N_reg] u;
  
  real<lower = 0> sigma_re;

  vector<lower = 0>[N_spl] sigma_B;
}

transformed parameters {
  vector[N_reg] convolved_re;
  
  vector[N] lp_tilde;
  vector[N] lp;
  
  vector[N] excessHaz;
  vector[N] cumExcessHaz;
  
  convolved_re = (sqrt(1 - rho) * v + sqrt(rho / scaling_factor) * u) * sigma_re;
  
  lp_tilde = linear_predictor(N, X_tilde, alpha);
  lp = linear_predictor_re(N, X, beta, region, convolved_re);
  
  excessHaz = hazLL(N, time .* exp(lp_tilde), mu, sigma, 0) .* exp(lp);
  cumExcessHaz = cumHazLL(N, time .* exp(lp_tilde), mu, sigma) .* exp(lp - lp_tilde);
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
      target += multi_normal_lpdf(beta[((i - 1) * df + 1):(i * df)] | rep_vector(0.0, df), pow(sigma_B[i], 2) * B[i]);
    }
  }

  // Linear Fixed coefficients
  for (i in 1:M_tilde) { 
    target += normal_lpdf(alpha[i] | 0, 10); 
  }
  
  for (i in (M_spl + 1):M) { 
    target += normal_lpdf(beta[i] | 0, 10); 
  }
  
  // LL location parameters
  target += normal_lpdf(mu | 0, 10);  
  
  // LL scale parameters
  target += cauchy_lpdf(sigma | 0, 1); 
  
  // Random effects
  target += icar_normal_1_lpdf(u | N_reg, node1, node2);
  target += normal_lpdf(v | 0, 1);
  
  // Hyperpriors
  target += normal_lpdf(sigma_re | 0, 1);
  target += beta_lpdf(rho | 0.5, 0.5);
  
  if (N_spl != 0) {
    for (i in 1:N_spl) {
      target += cauchy_lpdf(sigma_B[i] | 0, 1);
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
