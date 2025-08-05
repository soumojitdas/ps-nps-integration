// MRP model - traditional approach
// Models Y1 directly from demographics using NPS data
// PS data used only for estimating population cell counts N_c

data {
  // Sample sizes
  int<lower=0> N_nps;         // Non-probability sample size (Y1 observed)
  int<lower=1> J;             // Number of groups/states
  int<lower=1> K;             // Number of covariates (including intercept)
  
  // Group indicators (1-indexed)
  array[N_nps] int<lower=1, upper=J> group_nps;
  
  // Design matrices
  matrix[N_nps, K] X_nps;     // Covariates for NPS
  
  // Observed outcomes
  array[N_nps] int<lower=0, upper=1> Y1_nps;  // Outcome in NPS
  
  // Cell-level data for post-stratification
  int<lower=1> C;                // Number of post-stratification cells
  array[C] int<lower=1> n_ps_c;  // Observed PS sample sizes in each cell
  vector<lower=0>[C] w_bar_c;    // Mean survey weight in each cell
  array[C] int<lower=1, upper=J> cell_to_group; // Maps cells to groups
  
  // Cell mapping for NPS observations
  array[N_nps] int<lower=1, upper=C> nps_cell_id;  // Cell ID for each NPS observation
  
  // Cell-level design matrix
  matrix[C, K] X_cell;        // Design matrix for cell-level predictions
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
  
  // Likelihood for NPS: direct modeling of Y1
  // Linear predictor for NPS
  vector[N_nps] eta_prev_nps;
  for (i in 1:N_nps) {
    eta_prev_nps[i] = X_nps[i] * beta + u[group_nps[i]];
  }
  Y1_nps ~ bernoulli_logit(eta_prev_nps);
}

// // Generated quantities block for MRP WITH blending
// generated quantities {
//   // Calculate observed cell counts and means for NPS
//   array[C] int n_c_nps;
//   vector[C] y_bar_c_nps;
// 
//   // Initialize arrays
//   for (c in 1:C) {
//     n_c_nps[c] = 0;
//     y_bar_c_nps[c] = 0;
//   }
// 
//   // Count NPS observations and sums by cell
//   for (i in 1:N_nps) {
//     int c = nps_cell_id[i];
//     n_c_nps[c] += 1;
//     y_bar_c_nps[c] += Y1_nps[i];
//   }
// 
//   // Calculate cell means for NPS - handle cells with no observations
//   for (c in 1:C) {
//     if (n_c_nps[c] > 0) {
//       y_bar_c_nps[c] = y_bar_c_nps[c] / n_c_nps[c];
//     } else {
//       // If no observations, we'll use model prediction
//       y_bar_c_nps[c] = cell_prevalence[c];
//     }
//   }
// 
//   // Group-level MRP estimates using blending approach
//   vector[J] group_mrp_estimate;
//   real overall_mrp_estimate;
// 
//   // Calculate group-level MRP estimates with blending
//   for (g in 1:J) {
//     real total_estimate = 0;
//     real total_weight = 0;
// 
//     for (c in 1:C) {
//       if (cell_to_group[c] == g) {
//         // Blended estimate: use observed mean for sampled cells, model for unsampled
//         real cell_estimate;
// 
//         if (n_c_nps[c] > 0) {
//           // Use blending: observed for sampled units, model for unsampled
//           cell_estimate = n_c_nps[c] * y_bar_c_nps[c] +
//                          (N_c[c] - n_c_nps[c]) * cell_prevalence[c];
//         } else {
//           // No observations, use model prediction only
//           cell_estimate = N_c[c] * cell_prevalence[c];
//         }
// 
//         total_estimate += cell_estimate;
//         total_weight += N_c[c];
//       }
//     }
// 
//     // Only produce estimate if we have cells for this group
//     if (total_weight > 0) {
//       group_mrp_estimate[g] = total_estimate / total_weight;
//     } else {
//       group_mrp_estimate[g] = -1;  // Indicator for no data
//     }
//   }
// 
//   // Overall MRP estimate with blending
//   {
//     real total_estimate = 0;
//     real total_weight = 0;
// 
//     for (c in 1:C) {
//       // Blended cell estimate
//       real cell_estimate;
// 
//       if (n_c_nps[c] > 0) {
//         cell_estimate = n_c_nps[c] * y_bar_c_nps[c] +
//                        (N_c[c] - n_c_nps[c]) * cell_prevalence[c];
//       } else {
//         cell_estimate = N_c[c] * cell_prevalence[c];
//       }
// 
//       total_estimate += cell_estimate;
//       total_weight += N_c[c];
//     }
// 
//     if (total_weight > 0) {
//       overall_mrp_estimate = total_estimate / total_weight;
//     } else {
//       overall_mrp_estimate = -1;  // Indicator for no data
//     }
//   }
// }

// Generated quantities block for MRP WITHOUT blending
// This uses pure model predictions - standard MRP approach

generated quantities {
  // Cell-level predictions (model-based only)
  vector[C] cell_prev_Y1_model;

  // Calculate model predictions for each cell
  for (c in 1:C) {
    // Predict Y1 using the model (no direct observations)
    real logit_Y1 = X_cell[c] * beta + u[cell_to_group[c]];
    cell_prev_Y1_model[c] = inv_logit(logit_Y1);
  }

  // Group-level estimates using MODEL ONLY (no blending)
  vector[J] group_mrp_estimate;
  real overall_mrp_estimate;

  // Calculate group-level estimates
  for (g in 1:J) {
    real total_estimate = 0;
    real total_weight = 0;

    for (c in 1:C) {
      if (cell_to_group[c] == g) {
        // Use model prediction only (standard MRP)
        total_estimate += N_c[c] * cell_prev_Y1_model[c];
        total_weight += N_c[c];
      }
    }

    group_mrp_estimate[g] = total_weight > 0 ? total_estimate / total_weight : -1;
  }

  // Overall estimate (model-based only)
  {
    real total_estimate = 0;
    real total_weight = 0;

    for (c in 1:C) {
      // Pure model prediction
      total_estimate += N_c[c] * cell_prev_Y1_model[c];
      total_weight += N_c[c];
    }

    overall_mrp_estimate = total_estimate / total_weight;
  }
  
    // ===== POSTERIOR PREDICTIVE CHECKS =====
  // Generate replicated Y1 values for NPS sample
  array[N_nps] int Y1_nps_rep;
  vector[N_nps] prob_Y1_nps;
  
  for (i in 1:N_nps) {
    prob_Y1_nps[i] = inv_logit(X_nps[i] * beta + u[group_nps[i]]);
    Y1_nps_rep[i] = bernoulli_rng(prob_Y1_nps[i]);
  }
  
}
