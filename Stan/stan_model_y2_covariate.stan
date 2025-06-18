// Stan model for Y2 as predictor (not measurement error)
// Uses blending approach for NPS-based MRP estimates
// Y2 is generated first, Y1 depends on demographics and Y2

data {
  // Sample sizes
  int<lower=0> N_ps;          // PS: observe Y2, missing Y1
  int<lower=0> N_nps;         // NPS: observe Y1, missing Y2  
  int<lower=0> N_overlap;     // Overlap: observe both
  int<lower=1> J;             // Number of states/PUMAs/small areas
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
  
  // // Population cell counts
  // vector<lower=0>[C] N_c;
}

parameters {
  // Y2 model parameters 
  vector[K] gamma;              // Demographics -> Y2
  vector[J] v_raw;              // State random effects for Y2
  real<lower=0> sigma_v;        // SD of Y2 random effects
  
  // Y1 model parameters  
  vector[K] beta;               // Demographics -> Y1
  real alpha;                   // Y2 -> Y1 (KEY PARAMETER)
  vector[J] u_raw;              // State random effects for Y1
  real<lower=0> sigma_u;        // SD of Y1 random effects
  
  // Population cell counts
  vector<lower=0>[C] N_c;       // number of units in population cell c
}

transformed parameters {
  // Non-centered parameterization for random effects
  vector[J] v = sigma_v * v_raw;  // Y2 random effects
  vector[J] u = sigma_u * u_raw;  // Y1 random effects
  
  // Linear predictors for Y2 (demographics only)
  vector[N_ps] eta_Y2_ps;
  vector[N_nps] eta_Y2_nps;
  vector[N_overlap] eta_Y2_overlap;
  
  // Linear predictors for Y1 (demographics + Y2)
  // Note: add alpha*Y2 in the model block where Y2 is observed
  vector[N_ps] eta_Y1_base_ps;
  vector[N_nps] eta_Y1_base_nps;
  vector[N_overlap] eta_Y1_base_overlap;
  
  // Calculate Y2 lin preds
  for (i in 1:N_ps) {
    eta_Y2_ps[i] = X_ps[i] * gamma + v[group_ps[i]];
  }
  
  for (i in 1:N_nps) {
    eta_Y2_nps[i] = X_nps[i] * gamma + v[group_nps[i]];
  }
  
  for (i in 1:N_overlap) {
    eta_Y2_overlap[i] = X_overlap[i] * gamma + v[group_overlap[i]];
  }
  
  // Calculate Y1 base lin preds (without Y2 effect)
  for (i in 1:N_ps) {
    eta_Y1_base_ps[i] = X_ps[i] * beta + u[group_ps[i]];
  }
  
  for (i in 1:N_nps) {
    eta_Y1_base_nps[i] = X_nps[i] * beta + u[group_nps[i]];
  }
  
  for (i in 1:N_overlap) {
    eta_Y1_base_overlap[i] = X_overlap[i] * beta + u[group_overlap[i]];
  }
}

model {
  // Priors
  gamma ~ normal(0, 2);
  beta ~ normal(0, 2);
  alpha ~ normal(0, 2);  // Prior for Y2->Y1 effect
  
  v_raw ~ std_normal();
  u_raw ~ std_normal();
  sigma_v ~ normal(0, 0.5);
  sigma_u ~ normal(0, 0.5);
  
  // No prior on N_c (Stan defaults to flat prior)
  
  // Population cell size model
  for (c in 1:C) {
    n_ps_c[c] ~ poisson(N_c[c] / w_bar_c[c]);
  }
  
  // ===== PS SAMPLE =====
  // Observe Y2, model it directly
  Y2_ps ~ bernoulli_logit(eta_Y2_ps);
  
  // ===== NPS SAMPLE =====  
  // Observe Y1, but Y2 is missing
  // Thus, we need to marginalize over Y2
  for (i in 1:N_nps) {
    real p_Y2 = inv_logit(eta_Y2_nps[i]);
    real log_lik_Y2_0 = log(1 - p_Y2) + bernoulli_logit_lpmf(Y1_nps[i] | eta_Y1_base_nps[i] + alpha * 0);
    real log_lik_Y2_1 = log(p_Y2) + bernoulli_logit_lpmf(Y1_nps[i] | eta_Y1_base_nps[i] + alpha * 1);
    target += log_sum_exp(log_lik_Y2_0, log_lik_Y2_1);
  }
  
  // ===== OVERLAP SAMPLE =====
  // Observe both Y1 and Y2
  Y2_overlap ~ bernoulli_logit(eta_Y2_overlap);
  
  // // Y1 given Y2
  // for (i in 1:N_overlap) {
  //   Y1_overlap[i] ~ bernoulli_logit(eta_Y1_base_overlap[i] + alpha * Y2_overlap[i]);
  // }
  
  // Vectorized for efficiency
  Y1_overlap ~ bernoulli_logit(eta_Y1_base_overlap + alpha * to_vector(Y2_overlap));
}

generated quantities {
  // Impute Y1 for PS sample (where we observe Y2)
  vector[N_ps] prob_Y1_given_Y2;
  array[N_ps] int Y1_ps_imputed;
  
  for (i in 1:N_ps) {
    // Direct calculation: P(Y1=1|Y2,X) = logit^(-1)(X*beta + alpha*Y2 + u)
    prob_Y1_given_Y2[i] = inv_logit(eta_Y1_base_ps[i] + alpha * Y2_ps[i]);
    Y1_ps_imputed[i] = bernoulli_rng(prob_Y1_given_Y2[i]);
  }
  
  // Impute Y2 for NPS sample (where we observe Y1)
  vector[N_nps] prob_Y2_given_Y1;
  array[N_nps] int Y2_nps_imputed;
  
  for (i in 1:N_nps) {
    // Using Bayes rule: P(Y2=1|Y1,X) \prop P(Y1|Y2=1,X) * P(Y2=1|X)
    real p_Y2_prior = inv_logit(eta_Y2_nps[i]);
    real p_Y1_given_Y2_0 = inv_logit(eta_Y1_base_nps[i] + alpha * 0);
    real p_Y1_given_Y2_1 = inv_logit(eta_Y1_base_nps[i] + alpha * 1);
    
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
  
  // Add overlap observations to NPS counts and sums
  for (i in 1:N_overlap) {
    int c = overlap_cell_id[i];
    n_c_nps[c] += 1;
    y_bar_c_nps[c] += Y1_overlap[i];
  }
  
  // Calculate cell means for NPS - handle cells with no observations
  for (c in 1:C) {
    if (n_c_nps[c] > 0) {
      y_bar_c_nps[c] = y_bar_c_nps[c] / n_c_nps[c];
    } else {
      // If no observations, we'll use model prediction; won't happen in simulations
      y_bar_c_nps[c] = 0;  // Placeholder, will use model
    }
  }
  
  // Calculate cell-level predictions for cells without observations
  vector[C] cell_prev_Y2;
  vector[C] cell_prev_Y1_model;
  
  for (c in 1:C) {
    // Y2 prevalence in cell (for reference/diagnostics)
    cell_prev_Y2[c] = inv_logit(X_cell[c] * gamma + v[cell_to_group[c]]);
    
    // Y1 model prediction (marginalizing over Y2)
    real p_Y1_given_Y2_0 = inv_logit(X_cell[c] * beta + alpha * 0 + u[cell_to_group[c]]);
    real p_Y1_given_Y2_1 = inv_logit(X_cell[c] * beta + alpha * 1 + u[cell_to_group[c]]);
    cell_prev_Y1_model[c] = p_Y1_given_Y2_0 * (1 - cell_prev_Y2[c]) + 
                            p_Y1_given_Y2_1 * cell_prev_Y2[c];
  }
  
  // Group-level estimates using blending approach
  vector[J] group_mrp_estimate;
  real overall_mrp_estimate;
  
  // Calculate group-level MRP estimates
  for (g in 1:J) {
    real total_estimate = 0;
    real total_weight = 0;
    
    for (c in 1:C) {
      if (cell_to_group[c] == g) {
        // Blended estimate: use observed mean for sampled cells, model for unsampled
        real cell_estimate;
        
        if (n_c_nps[c] > 0) {
          // Use blending: observed for sampled units, model for unsampled
          cell_estimate = n_c_nps[c] * y_bar_c_nps[c] + 
                         (N_c[c] - n_c_nps[c]) * cell_prev_Y1_model[c];
        } else {
          // No observations, use model prediction only
          cell_estimate = N_c[c] * cell_prev_Y1_model[c];
        }
        
        total_estimate += cell_estimate;
        total_weight += N_c[c];
      }
    }
    
    group_mrp_estimate[g] = total_weight > 0 ? total_estimate / total_weight : -1;
  }
  
  // Overall MRP estimate
  {
    real total_estimate = 0;
    real total_weight = 0;

    for (c in 1:C) {
      // Blended cell estimate
      real cell_estimate;
      
      if (n_c_nps[c] > 0) {
        cell_estimate = n_c_nps[c] * y_bar_c_nps[c] + 
                       (N_c[c] - n_c_nps[c]) * cell_prev_Y1_model[c];
      } else {
        cell_estimate = N_c[c] * cell_prev_Y1_model[c];
      }
      
      total_estimate += cell_estimate;
      total_weight += N_c[c];
    }
    
    overall_mrp_estimate = total_estimate / total_weight;
  }
}
