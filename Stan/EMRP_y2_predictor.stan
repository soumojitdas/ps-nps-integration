// EMRP Stan model for Y2 as predictor (not measurement error)
// Age-only model matching the predictive setup
// Y2 is generated first, Y1 depends on age and Y2

data {
  // Sample sizes
  int<lower=0> N_ps;          // PS: observe Y2, missing Y1
  int<lower=0> N_nps;         // NPS: observe Y1, missing Y2  
  int<lower=0> N_overlap;     // Overlap: observe both
  int<lower=1> J;             // Number of states/PUMAs
  int<lower=1> K;             // Number of covariates (intercept + age categories)
  
  // Group indicators
  array[N_ps] int<lower=1, upper=J> group_ps;
  array[N_nps] int<lower=1, upper=J> group_nps;
  array[N_overlap] int<lower=1, upper=J> group_overlap;
  
  // Design matrices (age only)
  matrix[N_ps, K] X_ps;
  matrix[N_nps, K] X_nps;
  matrix[N_overlap, K] X_overlap;
  
  // Observed outcomes
  array[N_ps] int<lower=0, upper=1> Y2_ps;         // Y2 in PS
  array[N_nps] int<lower=0, upper=1> Y1_nps;       // Y1 in NPS
  array[N_overlap] int<lower=0, upper=1> Y1_overlap;
  array[N_overlap] int<lower=0, upper=1> Y2_overlap;
  
  // Cell-level data for population estimates
  int<lower=1> C;                // Number of cells (ST × age × Y2)
  array[C] int<lower=1> n_ps_c;  // Observed PS counts per cell
  vector<lower=0>[C] w_bar_c;    // Mean survey weight per cell
  array[C] int<lower=1, upper=J> cell_to_group; // Maps cells to states
  matrix[C, K] X_cell;           // Design matrix for cells
}

parameters {
  // Y2 model parameters 
  vector[K] gamma;              // Age effects -> Y2
  vector[J] v_raw;              // State random effects for Y2
  real<lower=0> sigma_v;        // SD of Y2 random effects
  
  // Y1 model parameters  
  vector[K] beta;               // Age effects -> Y1
  real alpha;                   // Y2 -> Y1 (KEY PARAMETER)
  vector[J] u_raw;              // State random effects for Y1
  real<lower=0> sigma_u;        // SD of Y1 random effects
  
  // Population cell counts
  vector<lower=0>[C] N_c;       // Population size in each cell
}

transformed parameters {
  // Non-centered parameterization for random effects
  vector[J] v = sigma_v * v_raw;  // Y2 random effects
  vector[J] u = sigma_u * u_raw;  // Y1 random effects
  
  // Linear predictors for Y2 (age only)
  vector[N_ps] eta_Y2_ps;
  vector[N_nps] eta_Y2_nps;
  vector[N_overlap] eta_Y2_overlap;
  
  // Linear predictors for Y1 (age only, Y2 effect added in model block)
  vector[N_ps] eta_Y1_base_ps;
  vector[N_nps] eta_Y1_base_nps;
  vector[N_overlap] eta_Y1_base_overlap;
  
  // Calculate Y2 linear predictors
  for (i in 1:N_ps) {
    eta_Y2_ps[i] = X_ps[i] * gamma + v[group_ps[i]];
  }
  
  for (i in 1:N_nps) {
    eta_Y2_nps[i] = X_nps[i] * gamma + v[group_nps[i]];
  }
  
  for (i in 1:N_overlap) {
    eta_Y2_overlap[i] = X_overlap[i] * gamma + v[group_overlap[i]];
  }
  
  // Calculate Y1 base linear predictors (without Y2 effect)
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
  gamma ~ normal(0, 2);         // Y2 age effects
  beta ~ normal(0, 2);          // Y1 age effects
  alpha ~ normal(0, 2);         // Y2->Y1 effect
  
  v_raw ~ std_normal();         // Y2 state effects
  u_raw ~ std_normal();         // Y1 state effects
  sigma_v ~ normal(0, 0.5);     // Y2 state SD
  sigma_u ~ normal(0, 0.5);     // Y1 state SD
  
  // Population cell size model
  for (c in 1:C) {
    n_ps_c[c] ~ poisson(N_c[c] / w_bar_c[c]);
  }
  
  // ===== PS SAMPLE =====
  // Observe Y2, model it directly
  Y2_ps ~ bernoulli_logit(eta_Y2_ps);
  
  // ===== NPS SAMPLE =====  
  // Observe Y1, but Y2 is missing
  // Marginalize over Y2
  for (i in 1:N_nps) {
    real p_Y2 = inv_logit(eta_Y2_nps[i]);
    real log_lik_Y2_0 = log(1 - p_Y2) + bernoulli_logit_lpmf(Y1_nps[i] | eta_Y1_base_nps[i] + alpha * 0);
    real log_lik_Y2_1 = log(p_Y2) + bernoulli_logit_lpmf(Y1_nps[i] | eta_Y1_base_nps[i] + alpha * 1);
    target += log_sum_exp(log_lik_Y2_0, log_lik_Y2_1);
  }
  
  // ===== OVERLAP SAMPLE =====
  // Observe both Y1 and Y2
  Y2_overlap ~ bernoulli_logit(eta_Y2_overlap);
  
  // Y1 given Y2 (vectorized)
  Y1_overlap ~ bernoulli_logit(eta_Y1_base_overlap + alpha * to_vector(Y2_overlap));
}

generated quantities {
  // EMRP: Only impute Y2 for NPS
  vector[N_nps] prob_Y2_given_Y1;
  array[N_nps] int Y2_nps_imputed;
  
  for (i in 1:N_nps) {
    // Using Bayes rule: P(Y2=1|Y1,X) ∝ P(Y1|Y2=1,X) * P(Y2=1|X)
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
  
  // That's it! Everything else happens in R:
  // - Cell assignment for NPS using imputed Y2
  // - Blended estimation
  // - Aggregation to state/overall levels
}
