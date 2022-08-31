functions {
  // -------------------------
  // Power Generalized Weibull
  // -------------------------

  // Hazard Function
  vector hazPGW (int N, vector time, real eta, real nu, real theta, int L) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = log(nu) - log(theta) - nu * log(eta) + (nu - 1) * log(time[i]) + ((1 / theta) - 1) * log(1 + pow((time[i] / eta), nu));
    }
    
    if (L == 1) {
      return res;
    } else {
      return exp(res);
    }
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

  // Liner Predictor
  vector linear_predictor (int N, matrix design_matrix, vector fixed_coeff, int[] region, vector random_effect) { 
    vector[N] res;
    vector[N] re_vector;
    for (i in 1:N) {
      re_vector[i] = random_effect[region[i]];
    }
    res = design_matrix * fixed_coeff + re_vector;
    
    return res;
  }
  
  // Distribution function for the ICAR model with constraint to the sum of u's
  // Reference: Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in Stan
  real icar_normal_lpdf (vector random_effect, int N_reg, int[] node1, int[] node2) { // For a general random effect "random_effect"
    return -0.5 * dot_self(random_effect[node1] - random_effect[node2]) + normal_lpdf(sum(random_effect) | 0, 0.001 * N_reg);
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
  
  vector[N_reg] u_tilde;
  vector[N_reg] u;
}

transformed parameters {
  vector[N] lp_tilde;
  vector[N] lp;
  
  real<lower=0> eta;
  real<lower=0> nu; 
  real<lower=0> theta;
  
  eta = exp(log_eta);
  nu = exp(log_nu);
  theta = exp(log_theta);
  
  lp_tilde = linear_predictor(N, X_tilde, alpha, region, u_tilde);
  lp = linear_predictor(N, X, beta, region, u);
}

model {
  // --------------
  // Log-likelihood
  // --------------
  
  {
    vector[N] excessHaz;
    vector[N] cumExcessHaz;
    
    excessHaz = hazPGW(N, time .* exp(lp_tilde), eta, nu, theta, 0) .* exp(lp);
    cumExcessHaz = cumHazPGW(N, time .* exp(lp_tilde), eta, nu, theta) .* exp(lp - lp_tilde);
    
    target += sum(log(pop_haz[cens] + excessHaz[cens])) - sum(cumExcessHaz);
  }
  
  // -------------------
  // Prior distributions
  // -------------------
  
  // Fixed coefficients
  alpha ~ normal(0, 10);
  beta  ~ normal(0, 10);
  
  // PGW scale parameters
  target += cauchy_lpdf(log_eta | 0, 2.5); 
  
  // PGW shape parameters
  target += cauchy_lpdf(log_nu | 0, 2.5);
  target += cauchy_lpdf(log_theta | 0, 2.5); // Check all the priors
  
  // Random effects
  target += icar_normal_lpdf(u_tilde | N_reg, node1, node2);
  target += icar_normal_lpdf(u | N_reg, node1, node2);
  
}

// generated quantities { } 
