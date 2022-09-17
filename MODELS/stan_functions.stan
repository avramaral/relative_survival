functions {
  // -------------------------------
  // Power Generalized Weibull (PGW)
  // -------------------------------

  // Hazard Function PGW
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
  
  // Cumulative Hazard Function PGW
  vector cumHazPGW (int N, vector time, real eta, real nu, real theta) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = - 1 + pow(1 + pow((time[i] / eta), nu), (1 / theta));
    }
    
    return res;
  } 
  
  // ---------------
  // Log Normal (LN)
  // ---------------

  // Hazard Function LN
  vector hazLN (int N, vector time, real mu, real sigma, int L) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = lognormal_lpdf(time[i] | mu, sigma) - lognormal_lccdf(time[i] | mu, sigma);
    }
    
    if (L == 1) {
      return res;
    } else {
      return exp(res);
    }
  }
  
  // Cumulative Hazard Function LN
  vector cumHazLN (int N, vector time, real mu, real sigma) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = - lognormal_lccdf(time[i] | mu, sigma);
    }

    return res;
  } 
  
  // -----------------
  // Log Logistic (LL)
  // -----------------

  // Hazard Function LL
  vector hazLL (int N, vector time, real mu, real sigma, int L) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = logistic_lpdf(log(time[i]) | mu, sigma) - log(time[i]) - logistic_lccdf(log(time[i]) | mu, sigma);
    }
    
    if (L == 1) {
      return res;
    } else {
      return exp(res);
    }
  }
  
  // Cumulative Hazard Function LL
  vector cumHazLL (int N, vector time, real mu, real sigma) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = - logistic_lccdf(log(time[i]) | mu, sigma);
    }
    
    return res;
  }
  
  // ----------------------
  // Generalized Gamma (GG)
  // ----------------------

  // Hazard Function GG
  vector hazGG (int N, vector time, real eta, real nu, real theta, int L) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = log(theta) - nu * log(eta) - lgamma(nu / theta) + (nu - 1) * log(time[i]) - pow(time[i] / eta, theta) - gamma_lccdf(pow(time[i], theta) | nu / theta, pow(eta, - theta));
    }
    
    if (L == 1) {
      return res;
    } else {
      return exp(res);
    }
  }
  
  // Cumulative Hazard Function GG
  vector cumHazGG (int N, vector time, real eta, real nu, real theta) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = - gamma_lccdf(pow(time[i], theta) | nu / theta, pow(eta, - theta));
    }
    
    return res;
  } 

  // -----------
  // Gamma (GAM)
  // -----------

  // Hazard Function GAM
  vector hazGAM (int N, vector time, real eta, real nu, int L) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = gamma_lpdf(time[i] | nu, 1 / eta) - gamma_lccdf(time[i] | nu, 1 / eta);
    }
    
    if (L == 1) {
      return res;
    } else {
      return exp(res);
    }
  }
  
  // Cumulative Hazard Function GAM
  vector cumHazGAM (int N, vector time, real eta, real nu) { 
    vector[N] res;
    for (i in 1:N) {
      res[i] = - gamma_lccdf(time[i] | nu, 1 / eta);
    }
    
    return res;
  } 
  
  // ------
  // Others
  // ------

  // Liner Predictor with no random effect
  vector linear_predictor (int N, matrix design_matrix, vector fixed_coeff) { 
    vector[N] res;
    res = design_matrix * fixed_coeff;
    
    return res;
  }

  // Liner Predictor with random effect
  vector linear_predictor_re (int N, matrix design_matrix, vector fixed_coeff, int[] region, vector random_effect) { 
    vector[N] res;
    vector[N] re_vector;
    for (i in 1:N) {
      re_vector[i] = random_effect[region[i]];
    }
    res = design_matrix * fixed_coeff + re_vector;
    
    return res;
  }
  
  // Distribution function for the ICAR model with constraint for the sum of u's = 0
  // Reference: Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in Stan
  real icar_normal_lpdf (vector random_effect, int N_reg, int[] node1, int[] node2) { // For a general random effect "random_effect"
    return -0.5 * dot_self(random_effect[node1] - random_effect[node2]) + normal_lpdf(sum(random_effect) | 0, 0.001 * N_reg);
  } 
  
}
