functions {
  // -------------------------
  // Power Generalized Weibull
  // -------------------------

  // Hazard Function
  vector hazPGW (int N, vector time, real eta, real nu, real theta) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = nu / (theta * pow(eta, nu)) * pow(time[i], (nu - 1)) * pow(1 + pow((time[i] / eta), nu), ((1 / theta) - 1));
    }
  
    return res;
  }
  
  // Cumulative Hazard Function
  vector cumHazPGW (int N, vector time, real eta, real nu, real theta) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = - 1 + pow(1 + pow((time[i] / eta), nu), (1 / theta));
    }
    
    return res;
  }

  // ------
  // Others
  // ------

  // Liner Predictor for X tilde
  vector linear_predictor_tilde (int N, matrix X_tilde, vector alpha) { // Need to include random effects.
    vector[N] res;
    res = X_tilde * alpha;
    
    return res;
  }
  
  // Linear Predictor for X
  vector linear_predictor (int N, matrix X, vector beta) { // Need to include random effects.
    vector[N] res;
    res = X * beta;
    
    return res;
  }
  
}

data {
  int<lower = 1> N;
  int<lower = 0> N_cens;
  int<lower = 0> M_tilde;
  int<lower = 0> M;
  int<lower = 1, upper = N> cens[N_cens];
  vector<lower = 0>[N] time;
  vector<lower = 0>[N] pop_haz;
  matrix[N, M_tilde] X_tilde;
  matrix[N, M] X;
  // Still missing information about the adjacency matrix.
}

parameters {
  vector[M_tilde] alpha;
  vector[M] beta;
  
  real<lower=0> eta;
  real<lower=0> nu;
  real<lower=0> theta; // Check the constraints.
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
    
    lp_tilde = linear_predictor_tilde (N, X_tilde, alpha);
    lp = linear_predictor (N, X, beta);
    
    excessHaz = hazPGW(N, time .* exp(lp_tilde), eta, nu, theta) .* exp(lp);
    cumExcessHaz = cumHazPGW(N, time .* exp(lp_tilde), eta, nu, theta) .* exp(lp - lp_tilde);
    
    target += sum(log(pop_haz[cens] + excessHaz[cens])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  // PGW scale parameters
  target += cauchy_lpdf(eta | 0, 1);
  
  // PGW shape parameters
  target += cauchy_lpdf(nu | 0, 1);
  target += gamma_lpdf(theta | 0.5, 0.5);
  
}
