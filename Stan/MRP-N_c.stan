data {
  // Sample sizes
  int<lower=0> N_ps;          // Probability sample size (Y2 observed, Y1 missing)
  int<lower=0> N_nps;         // Non-probability sample size (Y1 observed, but not used in this model)
  int<lower=1> J;             // Number of groups/states
  int<lower=1> K;             // Number of covariates (including intercept)
  
  // Group indicators (1-indexed)
  array[N_ps] int<lower=1, upper=J> group_ps;
  array[N_nps] int<lower=1, upper=J> group_nps;
  
  // Design matrices
  matrix[N_ps, K] X_ps;       // Covariates for PS
  matrix[N_nps, K] X_nps;     // Covariates for NPS
  
  // Observed outcomes
  array[N_nps] int<lower=0, upper=1> Y1_nps;  // Gold standard in NPS
  
  // Cell-level data for post-stratification
  int<lower=1> C;                // Number of post-stratification cells
  array[C] int<lower=1> n_ps_c;  // Observed PS sample sizes in each cell
  vector<lower=0>[C] w_bar_c;    // Mean survey weight in each cell
  array[C] int<lower=1, upper=J> cell_to_group; // Maps cells to groups
  
  // Cell-level design matrix
  matrix[C, K] X_cell;        // Design matrix for prevalence model
}

parameters {
  // Prevalence model parameters
  vector[K] beta;             // Regression coefficients for prevalence
  vector[J] u_raw;            // Raw group random effects for prevalence
  real<lower=0> sigma_u;      // SD of group random effects for prevalence
  
  // Population cell proportions 
  vector<lower=0>[C] N_c;     // Population counts in each cell
}

transformed parameters {
  // Group random effects (non-centered parameterization)
  vector[J] u = sigma_u * u_raw;  // Prevalence random effects
  
  // Cell-level predicted prevalence
  vector[C] eta_prev_cell;
  vector[C] cell_prevalence;
  
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
  sigma_u ~ cauchy(0, 2.5);
  
  // No prior on N_c - using flat prior
  
  // Poisson model for observed PS cell counts
  for (c in 1:C) {
    n_ps_c[c] ~ poisson(N_c[c] / w_bar_c[c]);
  }
  
  // Likelihood for NPS: direct modeling of gold standard (Y1)
  // Linear predictor for NPS
  vector[N_nps] eta_prev_nps;
  for (i in 1:N_nps) {
    eta_prev_nps[i] = X_nps[i] * beta + u[group_nps[i]];
  }
  Y1_nps ~ bernoulli_logit(eta_prev_nps);
}

generated quantities {
  // Group-level MRP estimates
  vector[J] group_mrp_estimate;
  real overall_mrp_estimate;
  
  // Calculate group-level MRP estimates
  for (g in 1:J) {
    real total_estimate = 0;
    real total_weight = 0;
    
    for (c in 1:C) {
      if (cell_to_group[c] == g) {
        // Use model prediction for each cell
        total_estimate += N_c[c] * cell_prevalence[c];
        total_weight += N_c[c];
      }
    }
    
    // Only produce estimate if we have cells for this group
    if (total_weight > 0) {
      group_mrp_estimate[g] = total_estimate / total_weight;
    } else {
      group_mrp_estimate[g] = -1;  // Indicator for no data
    }
  }
  
  // Overall MRP estimate
  {
    real total_estimate = 0;
    real total_weight = 0;

    for (c in 1:C) {
      total_estimate += N_c[c] * cell_prevalence[c];
      total_weight += N_c[c];
    }
    
    if (total_weight > 0) {
      overall_mrp_estimate = total_estimate / total_weight;
    } else {
      overall_mrp_estimate = -1;  // Indicator for no data
    }
  }
}
