// OG Model 3: Y1 = X*beta + (alpha + X*delta)*Y2 + gamma_j
// Global Y2 effect + demographic × Y2 interactions + marginal area effects
// Area effects are MARGINAL (not interacting with Y2)
// NO area random effects (u_j, v_j = 0)

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
  // Y2 model parameters
  vector[K] gamma;              // Demographics -> Y2
  
  // Y1 model parameters
  vector[K] beta;               // Demographics -> Y1
  real alpha;                   // Global Y2 -> Y1 effect
  vector[K-1] delta_raw;        // Demographic × Y2 interactions (excluding intercept)
  vector[J] gamma_j_raw;        // Marginal area effects (raw)
  real<lower=0> sigma_gamma;    // SD of marginal area effects
  
  // Population cell counts
  vector<lower=0>[C] N_c;       // number of units in population cell c
}

transformed parameters {
  // Marginal area effects (non-centered)
  vector[J] gamma_j = sigma_gamma * gamma_j_raw;
  
  // Demographic × Y2 interactions (delta[1] = 0 for intercept)
  vector[K] delta;
  delta[1] = 0;  // No interaction with intercept
  delta[2:K] = delta_raw;
  
  // Linear predictors for Y2 (demographics only, NO area effects)
  vector[N_ps] eta_Y2_ps;
  vector[N_nps] eta_Y2_nps;
  vector[N_overlap] eta_Y2_overlap;
  
  // Linear predictors for Y1 base (demographics + marginal area effects)
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
  
  // Calculate Y1 base linear predictors (including marginal area effects)
  for (i in 1:N_ps) {
    eta_Y1_base_ps[i] = X_ps[i] * beta + gamma_j[group_ps[i]];
  }
  
  for (i in 1:N_nps) {
    eta_Y1_base_nps[i] = X_nps[i] * beta + gamma_j[group_nps[i]];
  }
  
  for (i in 1:N_overlap) {
    eta_Y1_base_overlap[i] = X_overlap[i] * beta + gamma_j[group_overlap[i]];
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
    real p_Y2 = inv_logit(eta_Y2_nps[i]);
    
    // P(Y1=1|Y2=0) - includes marginal area effect
    real p_y1_given_y2_0 = inv_logit(eta_Y1_base_nps[i]);
    
    // P(Y1=1|Y2=1) - includes marginal area effect + Y2 effects
    real p_y1_given_y2_1 = inv_logit(eta_Y1_base_nps[i] + alpha + X_nps[i] * delta);
    
    real log_lik_Y2_0 = log(1 - p_Y2) + bernoulli_lpmf(Y1_nps[i] | p_y1_given_y2_0);
    real log_lik_Y2_1 = log(p_Y2) + bernoulli_lpmf(Y1_nps[i] | p_y1_given_y2_1);
    
    target += log_sum_exp(log_lik_Y2_0, log_lik_Y2_1);
  }
  
  // ===== OVERLAP SAMPLE =====
  // Observe both Y1 and Y2
  Y2_overlap ~ bernoulli_logit(eta_Y2_overlap);
  
  // Y1 given Y2 with marginal area effects and demographic interactions
  for (i in 1:N_overlap) {
    real eta_y1 = eta_Y1_base_overlap[i];  // Already includes marginal area effect
    
    if (Y2_overlap[i] == 1) {
      eta_y1 += alpha + X_overlap[i] * delta;
    }
    
    Y1_overlap[i] ~ bernoulli_logit(eta_y1);
  }
}

generated quantities {
  // Impute Y1 for PS sample (where we observe Y2)
  vector[N_ps] prob_Y1_given_Y2;
  array[N_ps] int Y1_ps_imputed;
  
  for (i in 1:N_ps) {
    real eta_y1 = eta_Y1_base_ps[i];  // Already includes marginal area effect
    
    if (Y2_ps[i] == 1) {
      eta_y1 += alpha + X_ps[i] * delta;
    }
    
    // Direct calculation: P(Y1=1|Y2,X)
    prob_Y1_given_Y2[i] = inv_logit(eta_y1);
    Y1_ps_imputed[i] = bernoulli_rng(prob_Y1_given_Y2[i]);
  }
  
  // Impute Y2 for NPS sample (where we observe Y1)
  vector[N_nps] prob_Y2_given_Y1;
  array[N_nps] int Y2_nps_imputed;
  
  for (i in 1:N_nps) {
    real p_Y2_prior = inv_logit(eta_Y2_nps[i]);
    real p_y1_given_y2_0 = inv_logit(eta_Y1_base_nps[i]);
    real p_y1_given_y2_1 = inv_logit(eta_Y1_base_nps[i] + alpha + X_nps[i] * delta);
    
    real log_posterior_Y2_0 = log(1 - p_Y2_prior) + bernoulli_lpmf(Y1_nps[i] | p_y1_given_y2_0);
    real log_posterior_Y2_1 = log(p_Y2_prior) + bernoulli_lpmf(Y1_nps[i] | p_y1_given_y2_1);
    
    prob_Y2_given_Y1[i] = exp(log_posterior_Y2_1 - log_sum_exp(log_posterior_Y2_0, log_posterior_Y2_1));
    Y2_nps_imputed[i] = bernoulli_rng(prob_Y2_given_Y1[i]);
  }
  
  // Impute Y1 for all cells for population estimates
  array[C] int Y1_cell_imputed;
  vector[C] prob_Y1_cell;
  vector[C] prob_Y2_cell;
  
  for (c in 1:C) {
    int j = cell_to_group[c];
    prob_Y2_cell[c] = inv_logit(X_cell[c] * gamma);
    
    // Marginal P(Y1=1) = P(Y1=1|Y2=0)*P(Y2=0) + P(Y1=1|Y2=1)*P(Y2=1)
    // Both include the marginal area effect gamma_j
    real p_y1_given_y2_0 = inv_logit(X_cell[c] * beta + gamma_j[j]);
    real p_y1_given_y2_1 = inv_logit(X_cell[c] * beta + gamma_j[j] + alpha + X_cell[c] * delta);
    
    prob_Y1_cell[c] = p_y1_given_y2_0 * (1 - prob_Y2_cell[c]) + 
                      p_y1_given_y2_1 * prob_Y2_cell[c];
    
    Y1_cell_imputed[c] = bernoulli_rng(prob_Y1_cell[c]);
  }
  
  // Group-level estimates (MRP-style)
  vector[J] group_mrp_estimate;
  real overall_mrp_estimate;
  
  for (j in 1:J) {
    real weighted_sum = 0;
    real total_weight = 0;
    
    for (c in 1:C) {
      if (cell_to_group[c] == j) {
        weighted_sum += N_c[c] * prob_Y1_cell[c];
        total_weight += N_c[c];
      }
    }
    
    group_mrp_estimate[j] = weighted_sum / total_weight;
  }
  
  // Overall estimate
  real weighted_sum_overall = 0;
  real total_weight_overall = 0;
  
  for (c in 1:C) {
    weighted_sum_overall += N_c[c] * prob_Y1_cell[c];
    total_weight_overall += N_c[c];
  }
  
  overall_mrp_estimate = weighted_sum_overall / total_weight_overall;
  
  // Store average total Y2 effect by area for diagnostics
  vector[J] avg_total_y2_effect_by_area;
  
  for (j in 1:J) {
    real weighted_sum_effect = 0;
    real total_weight_area = 0;
    
    for (c in 1:C) {
      if (cell_to_group[c] == j) {
        // Total Y2 effect for this cell: alpha + X*delta
        real cell_y2_effect = alpha + X_cell[c] * delta;
        weighted_sum_effect += N_c[c] * cell_y2_effect;
        total_weight_area += N_c[c];
      }
    }
    
    avg_total_y2_effect_by_area[j] = weighted_sum_effect / total_weight_area;
  }
}
