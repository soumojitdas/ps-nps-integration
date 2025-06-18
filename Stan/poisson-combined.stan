data {
  // Sample sizes
  int<lower=0> N_ps;          // Probability sample size (Y2 observed, Y1 missing)
  int<lower=0> N_nps;         // Non-probability sample size (Y1 observed, Y2 missing)
  int<lower=0> N_overlap;     // Overlap sample size (both Y1 and Y2 observed)
  int<lower=1> J;             // Number of groups/states
  int<lower=1> K;             // Number of covariates (including intercept)
  
  // Group indicators (1-indexed)
  array[N_ps] int<lower=1, upper=J> group_ps;
  array[N_nps] int<lower=1, upper=J> group_nps;
  array[N_overlap] int<lower=1, upper=J> group_overlap;
  
  // Design matrices
  matrix[N_ps, K] X_ps;          // Covariates for PS
  matrix[N_nps, K] X_nps;        // Covariates for NPS
  matrix[N_overlap, K] X_overlap; // Covariates for overlap
  
  // Observed outcomes
  array[N_ps] int<lower=0, upper=1> Y2_ps;         // Proxy in PS
  array[N_nps] int<lower=0, upper=1> Y1_nps;       // Gold standard in NPS
  array[N_overlap] int<lower=0, upper=1> Y1_overlap; // Gold standard in overlap
  array[N_overlap] int<lower=0, upper=1> Y2_overlap; // Proxy in overlap
  
  // Cell-level data for post-stratification
  int<lower=1> C;                // Number of post-stratification cells
  array[C] int<lower=1> n_ps_c;  // Observed PS sample sizes in each cell (from entire PS sample, including overlap)
  vector<lower=0>[C] w_bar_c;    // Mean survey weight in each cell
  array[C] int<lower=1, upper=J> cell_to_group; // Maps cells to groups
  
  // Cell mapping for observations
  array[N_nps] int<lower=1, upper=C> nps_cell_id;         // Cell ID for each NPS observation
  array[N_overlap] int<lower=1, upper=C> overlap_cell_id;  // Cell ID for each overlap observation
  
  // Cell-level design matrix
  matrix[C, K] X_cell;              // Design matrix for prevalence model
  
  // // Population cell proportions 
  // vector<lower=0>[C] N_c;
}

parameters {
  // Prevalence model parameters
  vector[K] beta;               // Regression coefficients for prevalence
  vector[J] u_raw;              // Raw group random effects for prevalence
  real<lower=0> sigma_u;        // SD of group random effects for prevalence
  
  // Measurement error parameters
  real alpha_sens;              // Sensitivity parameter (logit scale)
  real alpha_spec;              // Specificity parameter (logit scale)
  
  // Population cell proportions
  vector<lower=0>[C] N_c;
}

transformed parameters {
  // Group random effects (non-centered parameterization)
  vector[J] u = sigma_u * u_raw;  // Prevalence random effects

  // Sensitivity and specificity on probability scale
  real sens_prob = inv_logit(alpha_sens);
  real spec_prob = inv_logit(alpha_spec);
  
  // Linear predictors for prevalence
  vector[N_ps] eta_prev_ps;
  vector[N_nps] eta_prev_nps;
  vector[N_overlap] eta_prev_overlap;
  vector[C] eta_prev_cell;
  
  // Cell-level predicted prevalence
  vector[C] cell_prevalence;
  
  // Calculate prevalence linear predictors for individuals
  for (i in 1:N_ps) {
    eta_prev_ps[i] = X_ps[i] * beta + u[group_ps[i]];
  }
  
  for (i in 1:N_nps) {
    eta_prev_nps[i] = X_nps[i] * beta + u[group_nps[i]];
  }
  
  for (i in 1:N_overlap) {
    eta_prev_overlap[i] = X_overlap[i] * beta + u[group_overlap[i]];
  }
  
  // Calculate cell-level prevalence predictions
  for (c in 1:C) {
    eta_prev_cell[c] = X_cell[c] * beta + u[cell_to_group[c]];
    cell_prevalence[c] = inv_logit(eta_prev_cell[c]);
  }
}

model {
  // Priors for prevalence model
  beta ~ normal(0, 5);
  u_raw ~ std_normal();  // non-centered parameterization
  // sigma_u ~ normal(0, 1);
  sigma_u ~ cauchy(0, 2.5);
  // N_c ~ gamma(2, 0.002); // specifying nothing will induce flat-prior on N_c
  
  // Priors for measurement error parameters
  alpha_sens ~ normal(0, 1.5);
  alpha_spec ~ normal(0, 1.5);
  
  // Poisson model for observed PS cell counts from ENTIRE PS sample (including overlap)
  for (c in 1:C) {
    n_ps_c[c] ~ poisson(N_c[c] / w_bar_c[c]);
  }
  
  // Likelihood for NPS: direct modeling of gold standard (Y1)
  Y1_nps ~ bernoulli_logit(eta_prev_nps);
  
  // Likelihood for PS: modeling error-prone proxy (Y2) conditionally
  for (i in 1:N_ps) {
    real p_Y1 = inv_logit(eta_prev_ps[i]);
    
    if (Y2_ps[i] == 1) {
      // P(Y2=1) = sens*P(Y1=1) + (1-spec)*P(Y1=0)
      target += log_mix(p_Y1, 
                      bernoulli_lpmf(1 | sens_prob),
                      bernoulli_lpmf(1 | 1 - spec_prob));
    } else {
      // P(Y2=0) = (1-sens)*P(Y1=1) + spec*P(Y1=0)
      target += log_mix(p_Y1, 
                      bernoulli_lpmf(0 | sens_prob),
                      bernoulli_lpmf(0 | 1 - spec_prob));
    }
  }
  
  // Likelihood for overlap sample
  // Y1 follows prevalence model
  Y1_overlap ~ bernoulli_logit(eta_prev_overlap);
  
  // Y2 follows conditional model based on Y1
  for (i in 1:N_overlap) {
    if (Y1_overlap[i] == 1) {
      // For true positives, use sensitivity
      Y2_overlap[i] ~ bernoulli(sens_prob);
    } else {
      // For true negatives, use (inverse) specificity
      Y2_overlap[i] ~ bernoulli(1 - spec_prob);
    }
  }
}

generated quantities {
  // Impute gold standard Y1 for PS sample
  array[N_ps] int Y1_ps_imputed;
  vector[N_ps] prob_Y1_given_Y2;
  
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
  
  // Add overlap observations to counts and sums
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
      // If no observations, use model prediction
      y_bar_c_nps[c] = cell_prevalence[c];
    }
  }
  
  // Group-level MRP estimates using blended approach
  vector[J] group_mrp_estimate;
  real overall_mrp_estimate;
  
  // Calculate group-level MRP estimates
  for (g in 1:J) {
    real total_estimate = 0;
    real total_weight = 0;
    
    for (c in 1:C) {
      if (cell_to_group[c] == g) {
        // Calculate blended cell estimate
        real cell_estimate;
        
        if (n_c_nps[c] > 0) {
          cell_estimate = n_c_nps[c] * y_bar_c_nps[c] + 
                         (N_c[c] - n_c_nps[c]) * cell_prevalence[c];
        } else {
          cell_estimate = N_c[c] * cell_prevalence[c];
        }
        
        total_estimate += cell_estimate;
        total_weight += N_c[c];
      }
    }
    
    group_mrp_estimate[g] = total_estimate / total_weight;
  }
  
  // Overall MRP estimate
  {
    real total_estimate = 0;
    real total_weight = 0;

    for (c in 1:C) {
      // Calculate blended cell estimate
      real cell_estimate;
      
      if (n_c_nps[c] > 0) {
        cell_estimate = n_c_nps[c] * y_bar_c_nps[c] + 
                      (N_c[c] - n_c_nps[c]) * cell_prevalence[c];
      } else {
        cell_estimate = N_c[c] * cell_prevalence[c];
      }
      
      total_estimate += cell_estimate;
      total_weight += N_c[c];
    }
    
    overall_mrp_estimate = total_estimate / total_weight;
  }
  
  // Calculate P(Y1=1|Y2) for each PS individual
  for (i in 1:N_ps) {
    real p_Y1 = inv_logit(eta_prev_ps[i]);
    
    real numerator, denominator;
    
    if (Y2_ps[i] == 1) {
      // P(Y1=1|Y2=1) 
      numerator = sens_prob * p_Y1;
      denominator = sens_prob * p_Y1 + (1 - spec_prob) * (1 - p_Y1);
    } else {
      // P(Y1=1|Y2=0)
      numerator = (1 - sens_prob) * p_Y1;
      denominator = (1 - sens_prob) * p_Y1 + spec_prob * (1 - p_Y1);
    }
    
    prob_Y1_given_Y2[i] = numerator / denominator;
    
    // Sample Y1 given Y2
    Y1_ps_imputed[i] = bernoulli_rng(prob_Y1_given_Y2[i]);
  }
}
