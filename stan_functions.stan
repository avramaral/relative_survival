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
  
  // Distribution function for the ICAR model with constraint to the sum of u's
  // Reference: Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in Stan
  real icar_normal_lpdf (vector random_effect, int N_reg, int[] node1, int[] node2) { // For a general random effect "random_effect"
    return -0.5 * dot_self(random_effect[node1] - random_effect[node2]) + normal_lpdf(sum(random_effect) | 0, 0.001 * N_reg);
  } 
  
}
