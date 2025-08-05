// OG model following stan_model_y2_covariate.stan structure
// Enhanced with area-specific Y2 effects and demographic × Y2 interactions
// NO area random effects (u_j, v_j = 0)
// NOW STORES TOTAL Y2 EFFECTS FOR COMPARISON

data {
  // Sample sizes
  int<lower=0> N_ps;          // PS: observe Y2, missing Y1
  int<lower=0> N_nps;         // NPS: observe Y1, missing Y2  
  int<lower=0> N_overlap;     // Overlap: observe both
  int<lower=1> J;             // Number of states/areas
  int<lower=1> K;             // Number of demographic covariates
  
  // Group indicators
  array[N_ps] int<lower=1, upper=J> group_ps;
  array[N_nps] int<lower=1, upper=J> group_nps;
  array[N_overlap] int<lower=1, upper=J> group_overlap;
  
  // Design matrices (demographics only)
  matrix[N_ps, K] X_ps;
  matrix[N_nps, K] X_nps;
  matrix[N_overlap, K] X_overlap;
  
  // Observed outcomes
  array[N_ps] int<lower=0, upper=1> Y2_ps;         // Y2 in PS
  array[N_nps] int<lower=0, upper=1> Y1_nps;       // Y1 in NPS
  array[N_overlap] int<lower=0, upper=1> Y1_overlap;
  array[N_overlap] int<lower=0, upper=1> Y2_overlap;
  
  // Cell-level data for population estimates
  int<lower=1> C;
  array[C] int<lower=1> n_ps_c;
  vector<lower=0>[C] w_bar_c;
  array[C] int<lower=1, upper=J> cell_to_group;
  matrix[C, K] X_cell;
  
  // Cell mappings
  array[N_nps] int<lower=1, upper=C> nps_cell_id;
  array[N_overlap] int<lower=1, upper=C> overlap_cell_id;
}

parameters {
  // Y2 model parameters (NO v_j since sigma_v = 0)
  vector[K] gamma;              // Demographics -> Y2
  
  // Y1 model parameters (NO u_j since sigma_u = 0)
  vector[K] beta;               // Demographics -> Y1
  real alpha;                   // Global Y2 -> Y1 effect
  vector[K-1] delta_raw;        // Demographic × Y2 interactions (excluding intercept)
  vector[J] gamma_j_raw;        // Area-specific Y2 effects (raw)
  real<lower=0> sigma_gamma;    // SD of area-specific Y2 effects
  
  // Population cell counts
  vector<lower=0>[C] N_c;       // number of units in population cell c
}

transformed parameters {
  // Area-specific Y2 effects (non-centered)
  vector[J] gamma_j = sigma_gamma * gamma_j_raw;
  
  // Demographic × Y2 interactions (delta[1] = 0 for intercept)
  vector[K] delta;
  delta[1] = 0;  // No interaction with intercept
  delta[2:K] = delta_raw;
  
  // Base Y2 effect by area (alpha + gamma_j only)
  vector[J] base_y2_effect_by_area = alpha + gamma_j;
  
  // Linear predictors for Y2 (demographics only, NO area effects)
  vector[N_ps] eta_Y2_ps;
  vector[N_nps] eta_Y2_nps;
  vector[N_overlap] eta_Y2_overlap;
  
  // Linear predictors for Y1 base (demographics only, NO area effects)
  vector[N_ps] eta_Y1_base_ps;
  vector[N_nps] eta_Y1_base_nps;
  vector[N_overlap] eta_Y1_base_overlap;
  
  // Calculate Y2 linear predictors
  for (i in 1:N_ps) {
    eta_Y2_ps[i] = X_ps[i] * gamma;
  }
  
  for (i in 1:N_nps) {
    eta_Y2_nps[i] = X_nps[i] * gamma;
  }
  
  for (i in 1:N_overlap) {
    eta_Y2_overlap[i] = X_overlap[i] * gamma;
  }
  
  // Calculate Y1 base linear predictors
  for (i in 1:N_ps) {
    eta_Y1_base_ps[i] = X_ps[i] * beta;
  }
  
  for (i in 1:N_nps) {
    eta_Y1_base_nps[i] = X_nps[i] * beta;
  }
  
  for (i in 1:N_overlap) {
    eta_Y1_base_overlap[i] = X_overlap[i] * beta;
  }
}

model {
  // Priors
  gamma ~ normal(0, 2);
  beta ~ normal(0, 2);
  alpha ~ normal(0, 2);
  delta_raw ~ normal(0, 1);
  gamma_j_raw ~ std_normal();
  sigma_gamma ~ normal(0, 2);
  
  // No prior on N_c (Stan defaults to flat prior)
  
  // Population cell size model
  for (c in 1:C) {
    n_ps_c[c] ~ poisson(N_c[c] / w_bar_c[c]);
  }
  
  // ===== PS SAMPLE =====
  // Observe Y2, model it directly
  Y2_ps ~ bernoulli_logit(eta_Y2_ps);
  
  // ===== NPS SAMPLE =====  
  // Observe Y1, but Y2 is missing - marginalize over Y2
  for (i in 1:N_nps) {
    int j = group_nps[i];
    real p_Y2 = inv_logit(eta_Y2_nps[i]);
    
    // P(Y1=1|Y2=0)
    real p_y1_given_y2_0 = inv_logit(eta_Y1_base_nps[i]);
    
    // P(Y1=1|Y2=1) with area effect and interactions
    real p_y1_given_y2_1 = inv_logit(eta_Y1_base_nps[i] + 
                                      base_y2_effect_by_area[j] + 
                                      X_nps[i] * delta);
    
    real log_lik_Y2_0 = log(1 - p_Y2) + bernoulli_lpmf(Y1_nps[i] | p_y1_given_y2_0);
    real log_lik_Y2_1 = log(p_Y2) + bernoulli_lpmf(Y1_nps[i] | p_y1_given_y2_1);
    
    target += log_sum_exp(log_lik_Y2_0, log_lik_Y2_1);
  }
  
  // ===== OVERLAP SAMPLE =====
  // Observe both Y1 and Y2
  Y2_overlap ~ bernoulli_logit(eta_Y2_overlap);
  
  // Y1 given Y2 with area effects and interactions
  for (i in 1:N_overlap) {
    int j = group_overlap[i];
    real eta_y1 = eta_Y1_base_overlap[i];
    
    if (Y2_overlap[i] == 1) {
      eta_y1 += base_y2_effect_by_area[j] + X_overlap[i] * delta;
    }
    
    Y1_overlap[i] ~ bernoulli_logit(eta_y1);
  }
}

generated quantities {
  // Impute Y1 for PS sample (where we observe Y2)
  vector[N_ps] prob_Y1_given_Y2;
  array[N_ps] int Y1_ps_imputed;
  
  for (i in 1:N_ps) {
    int j = group_ps[i];
    real eta_y1 = eta_Y1_base_ps[i];
    
    if (Y2_ps[i] == 1) {
      eta_y1 += base_y2_effect_by_area[j] + X_ps[i] * delta;
    }
    
    // Direct calculation: P(Y1=1|Y2,X)
    prob_Y1_given_Y2[i] = inv_logit(eta_y1);
    Y1_ps_imputed[i] = bernoulli_rng(prob_Y1_given_Y2[i]);
  }
  
  // Impute Y2 for NPS sample (where we observe Y1)
  vector[N_nps] prob_Y2_given_Y1;
  array[N_nps] int Y2_nps_imputed;
  
  for (i in 1:N_nps) {
    int j = group_nps[i];
    // Using Bayes rule: P(Y2=1|Y1,X) ∝ P(Y1|Y2=1,X) * P(Y2=1|X)
    real p_Y2_prior = inv_logit(eta_Y2_nps[i]);
    real p_Y1_given_Y2_0 = inv_logit(eta_Y1_base_nps[i]);
    real p_Y1_given_Y2_1 = inv_logit(eta_Y1_base_nps[i] + 
                                      base_y2_effect_by_area[j] + 
                                      X_nps[i] * delta);
    
    real lik_Y2_0, lik_Y2_1;
    if (Y1_nps[i] == 1) {
      lik_Y2_0 = p_Y1_given_Y2_0 * (1 - p_Y2_prior);
      lik_Y2_1 = p_Y1_given_Y2_1 * p_Y2_prior;
    } else {
      lik_Y2_0 = (1 - p_Y1_given_Y2_0) * (1 - p_Y2_prior);
      lik_Y2_1 = (1 - p_Y1_given_Y2_1) * p_Y2_prior;
    }
    
    prob_Y2_given_Y1[i] = lik_Y2_1 / (lik_Y2_0 + lik_Y2_1);
    Y2_nps_imputed[i] = bernoulli_rng(prob_Y2_given_Y1[i]);
  }
  
  // Calculate observed cell counts and means for NPS
  array[C] int n_c_nps;
  vector[C] y_bar_c_nps;
  
  // Initialize arrays
  for (c in 1:C) {
    n_c_nps[c] = 0;
    y_bar_c_nps[c] = 0;
  }
  
  // Count NPS observations and sums by cell
  for (i in 1:N_nps) {
    int c = nps_cell_id[i];
    n_c_nps[c] += 1;
    y_bar_c_nps[c] += Y1_nps[i];
  }
  
  // Add overlap observations to NPS counts
  for (i in 1:N_overlap) {
    int c = overlap_cell_id[i];
    n_c_nps[c] += 1;
    y_bar_c_nps[c] += Y1_overlap[i];
  }
  
  // Calculate cell means
  for (c in 1:C) {
    if (n_c_nps[c] > 0) {
      y_bar_c_nps[c] = y_bar_c_nps[c] / n_c_nps[c];
    } else {
      // If no observations, use model prediction
      // Marginalize over Y2
      int j = cell_to_group[c];
      real p_y2 = inv_logit(X_cell[c] * gamma);
      real p_y1_given_y2_0 = inv_logit(X_cell[c] * beta);
      real p_y1_given_y2_1 = inv_logit(X_cell[c] * beta + 
                                        base_y2_effect_by_area[j] + 
                                        X_cell[c] * delta);
      y_bar_c_nps[c] = p_y1_given_y2_0 * (1 - p_y2) + p_y1_given_y2_1 * p_y2;
    }
  }
  
  // Group-level MRP estimates (blended approach)
  vector[J] group_mrp_estimate;
  real overall_mrp_estimate;
  
  // Calculate group estimates using blended approach
  for (j in 1:J) {
    real total_direct = 0;
    real total_model = 0;
    real total_weight = 0;
    
    for (c in 1:C) {
      if (cell_to_group[c] == j) {
        if (n_c_nps[c] > 0) {
          // Blended estimate
          total_direct += n_c_nps[c] * y_bar_c_nps[c];
          
          // Model prediction for unobserved units
          real p_y2 = inv_logit(X_cell[c] * gamma);
          real p_y1_given_y2_0 = inv_logit(X_cell[c] * beta);
          real p_y1_given_y2_1 = inv_logit(X_cell[c] * beta + 
                                            base_y2_effect_by_area[j] + 
                                            X_cell[c] * delta);
          real cell_prev = p_y1_given_y2_0 * (1 - p_y2) + p_y1_given_y2_1 * p_y2;
          
          total_model += (N_c[c] - n_c_nps[c]) * cell_prev;
          total_weight += N_c[c];
        } else {
          // Pure model estimate
          real p_y2 = inv_logit(X_cell[c] * gamma);
          real p_y1_given_y2_0 = inv_logit(X_cell[c] * beta);
          real p_y1_given_y2_1 = inv_logit(X_cell[c] * beta + 
                                            base_y2_effect_by_area[j] + 
                                            X_cell[c] * delta);
          real cell_prev = p_y1_given_y2_0 * (1 - p_y2) + p_y1_given_y2_1 * p_y2;
          
          total_model += N_c[c] * cell_prev;
          total_weight += N_c[c];
        }
      }
    }
    
    group_mrp_estimate[j] = total_weight > 0 ? (total_direct + total_model) / total_weight : -1;
  }
  
  // Overall estimate
  {
    real total_direct = 0;
    real total_model = 0; 
    real total_weight = 0;

    for (c in 1:C) {
      if (n_c_nps[c] > 0) {
        total_direct += n_c_nps[c] * y_bar_c_nps[c];
        
        // Model prediction for unobserved units
        int j = cell_to_group[c];
        real p_y2 = inv_logit(X_cell[c] * gamma);
        real p_y1_given_y2_0 = inv_logit(X_cell[c] * beta);
        real p_y1_given_y2_1 = inv_logit(X_cell[c] * beta + 
                                          base_y2_effect_by_area[j] + 
                                          X_cell[c] * delta);
        real cell_prev = p_y1_given_y2_0 * (1 - p_y2) + p_y1_given_y2_1 * p_y2;
        
        total_model += (N_c[c] - n_c_nps[c]) * cell_prev;
        total_weight += N_c[c];
      } else {
        // Pure model estimate
        int j = cell_to_group[c];
        real p_y2 = inv_logit(X_cell[c] * gamma);
        real p_y1_given_y2_0 = inv_logit(X_cell[c] * beta);
        real p_y1_given_y2_1 = inv_logit(X_cell[c] * beta + 
                                          base_y2_effect_by_area[j] + 
                                          X_cell[c] * delta);
        real cell_prev = p_y1_given_y2_0 * (1 - p_y2) + p_y1_given_y2_1 * p_y2;
        
        total_model += N_c[c] * cell_prev;
        total_weight += N_c[c];
      }
    }
    
    overall_mrp_estimate = (total_direct + total_model) / total_weight;
  }
  
  // ===== NEW: CALCULATE AVERAGE TOTAL Y2 EFFECTS BY AREA =====
  // This includes alpha + gamma_j + average(X*delta) for each area
  vector[J] avg_total_y2_effect_by_area;
  
  for (j in 1:J) {
    real sum_effect = 0;
    real total_pop = 0;  // Use real type for population count
    
    // Average over all cells in this area
    for (c in 1:C) {
      if (cell_to_group[c] == j) {
        // Complete effect for this cell: α + γ_j + X_cell*δ
        real cell_effect = alpha + gamma_j[j] + X_cell[c] * delta;
        sum_effect += cell_effect * N_c[c];  // Weight by population size
        total_pop += N_c[c];
      }
    }
    
    avg_total_y2_effect_by_area[j] = total_pop > 0 ? sum_effect / total_pop : 0;
  }
  
  // ===== ALSO STORE INDIVIDUAL PARAMETERS FOR DIAGNOSTICS =====
  real estimated_alpha = alpha;
  real estimated_sigma_gamma = sigma_gamma;
  vector[K] estimated_beta = beta;
  vector[K] estimated_gamma = gamma;
  vector[K] estimated_delta = delta;
  vector[J] estimated_gamma_j = gamma_j;
}
