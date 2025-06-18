library(cmdstanr)
library(rstanarm)
library(posterior)
library(tidyverse)
library(future)
library(future.apply)
library(survey)

# Helper function for survey design
create_survey_design <- function(data, sampling_scheme) {
  if(sampling_scheme == "Stratified") {
    svydesign(
      ids     = ~1,
      strata  = ~age_category,
      weights = ~survey_weight,
      data    = data
    )
  } else {
    # For SRS and PO, no strata
    svydesign(
      ids     = ~1,
      weights = ~survey_weight,
      data    = data
    )
  }
}

# RECOMMENDED: Randomized Realistic Selection
# ============================================
# Draw coefficients that ensure realistic patterns without extreme values

generate_realistic_nps_selection <- function(scenario = "mixed", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (scenario == "internet") {
    # Internet panel: monotonic decrease with age
    beta_0 <- rnorm(1, mean = -2.2, sd = 0.2)
    beta_age <- c(0, 
                  runif(1, -0.3, -0.1),   # 30-44
                  runif(1, -0.5, -0.2),   # 45-59
                  runif(1, -0.7, -0.3),   # 60-74
                  runif(1, -0.9, -0.5))   # 75+
    
  } else if (scenario == "health") {
    # Health survey: increase then slight decrease
    beta_0 <- rnorm(1, mean = -2.5, sd = 0.2)
    beta_age <- c(0,
                  runif(1, 0.0, 0.3),    # 30-44
                  runif(1, 0.2, 0.5),    # 45-59
                  runif(1, 0.3, 0.6),    # 60-74
                  runif(1, 0.1, 0.4))    # 75+
    
  } else if (scenario == "community") {
    # Community survey: peak at middle age
    beta_0 <- rnorm(1, mean = -2.3, sd = 0.2)
    beta_age <- c(0,
                  runif(1, 0.2, 0.5),    # 30-44
                  runif(1, 0.2, 0.5),    # 45-59
                  runif(1, -0.1, 0.2),   # 60-74
                  runif(1, -0.5, -0.1))  # 75+
    
  } else {  # "mixed" - most realistic
    # Mixed patterns: small random deviations
    beta_0 <- rnorm(1, mean = -2.3, sd = 0.15)
    beta_age <- c(0, rnorm(4, mean = 0, sd = 0.25))
    
    # Ensure no extreme values
    beta_age <- pmax(pmin(beta_age, 0.6), -0.6)
  }
  
  return(list(
    beta_0 = beta_0,
    beta_age = beta_age
  ))
}

# For a single simulation study:
nps_params <- generate_realistic_nps_selection(scenario = "health", seed = 17)
nps_beta_0 <- nps_params$beta_0
nps_beta_age <- nps_params$beta_age

create_overlap_samples_EMRP_predictive <- function(pop_data,
                                                   N_ps = 1100,
                                                   N_nps = 4300,
                                                   overlap_pct = 0.20,
                                                   sampling_scheme = "SRS",  # Options: "SRS", "Stratified", "PO"
                                                   # PS selection parameters (full vector with intercept first)
                                                   ps_betas = c(-3.1, 0.2, 0.3, 0.4, 0.2),  # length 5: intercept + 4 age effects
                                                   ps_sigma = 0.1,
                                                   # NPS selection parameters (full vector with intercept first)
                                                   nps_betas = c(-2.703002, 0.1404789, 0.4330459, 0.4223657, 0.2616391),  # length 5
                                                   nps_sigma = 0.1,
                                                   # MNAR mechanism parameters
                                                   use_mnar = FALSE,
                                                   mnar_shift = 0.5,
                                                   seed = 17) {
  set.seed(seed)
  
  # Validate inputs
  required_vars <- c("Y_1", "Y_2", "ST", "age_category")
  missing_vars <- setdiff(required_vars, names(pop_data))
  if (length(missing_vars) > 0) {
    stop("Population data missing required variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Validate beta vectors
  if (length(ps_betas) != 5) {
    stop("ps_betas must have length 5 (intercept + 4 age effects)")
  }
  if (length(nps_betas) != 5) {
    stop("nps_betas must have length 5 (intercept + 4 age effects)")
  }
  
  # Verify sampling scheme
  sampling_scheme <- match.arg(sampling_scheme, c("SRS", "Stratified", "PO"))
  
  # Ensure age factor is properly set
  pop_data <- pop_data %>%
    mutate(
      age_category = factor(age_category,
                            levels = c("18-29", "30-44", "45-59", "60-74", "75+"))
    )
  
  cat("=== EMRP PREDICTIVE SAMPLING ===\n")
  cat("Sampling scheme:", sampling_scheme, "\n")
  cat("MNAR mechanism:", ifelse(use_mnar, "ENABLED", "DISABLED"), "\n")
  cat("Using PS betas:", round(ps_betas, 3), "\n")
  cat("Using NPS betas:", round(nps_betas, 3), "\n")
  
  # =====================================================
  # PROBABILITY SAMPLE (PS) SELECTION
  # =====================================================
  
  if (sampling_scheme == "SRS") {
    # Simple Random Sampling
    ps_set <- sample(1:nrow(pop_data), N_ps, replace = FALSE)
    cat("\nPS size from SRS:", length(ps_set), "\n")
    
    total_pop_size <- nrow(pop_data)
    ps_weight_value <- total_pop_size / N_ps
    
    assign_ps_weights <- function(data) {
      data %>% mutate(survey_weight = ps_weight_value)
    }
    
  } else if (sampling_scheme == "Stratified") {
    # Stratified by age with proportional allocation
    age_strata <- pop_data %>%
      group_by(age_category) %>%
      summarise(N_h = n(), .groups = "drop") %>%
      mutate(n_h = round(N_ps * N_h / sum(N_h)))
    
    # Adjust to ensure exactly N_ps
    diff <- N_ps - sum(age_strata$n_h)
    if (diff != 0) age_strata$n_h[1] <- age_strata$n_h[1] + diff
    
    cat("\nStratified sampling allocation by age:\n")
    print(age_strata)
    
    ps_set <- integer(0)
    for (i in 1:nrow(age_strata)) {
      stratum <- age_strata$age_category[i]
      target_n <- age_strata$n_h[i]
      stratum_indices <- which(pop_data$age_category == stratum)
      
      if (target_n > length(stratum_indices)) {
        warning(sprintf("Requested %d samples from stratum %s, but only %d units available",
                        target_n, stratum, length(stratum_indices)))
        sampled <- stratum_indices
      } else {
        sampled <- sample(stratum_indices, target_n, replace = FALSE)
      }
      ps_set <- c(ps_set, sampled)
    }
    
    cat("Actual PS sample size:", length(ps_set), "(target:", N_ps, ")\n")
    
    assign_ps_weights <- function(data) {
      data %>%
        left_join(age_strata %>% select(age_category, N_h, n_h), by = "age_category") %>%
        mutate(survey_weight = N_h / n_h) %>%
        select(-N_h, -n_h)
    }
    
  } else if (sampling_scheme == "PO") {
    # Poisson sampling using age
    cat("\nPS Poisson sampling parameters:\n")
    cat("ps_betas:", round(ps_betas, 3), "\n")
    cat("ps_sigma:", ps_sigma, "\n")
    
    # Create design matrix WITH intercept
    X_ps_selection <- model.matrix(~ age_category, data = pop_data)
    
    # Validate dimensions
    if (ncol(X_ps_selection) != length(ps_betas)) {
      stop("Mismatch: design matrix has ", ncol(X_ps_selection), 
           " columns but ps_betas has length ", length(ps_betas))
    }
    
    # Generate errors and calculate selection probabilities
    ps_errors <- rnorm(nrow(pop_data), 0, ps_sigma)
    logit_p <- X_ps_selection %*% ps_betas + ps_errors
    selection_prob <- plogis(logit_p)
    
    # Report selection probabilities by age
    cat("\nPS selection probabilities by age:\n")
    ps_prob_summary <- pop_data %>%
      mutate(selection_prob = as.vector(selection_prob)) %>%
      group_by(age_category) %>%
      summarise(
        mean_prob = mean(selection_prob),
        sd_prob = sd(selection_prob),
        .groups = "drop"
      )
    print(ps_prob_summary)
    
    # Poisson sampling
    pop_data <- pop_data %>%
      mutate(
        ps_selection_prob = as.vector(selection_prob),
        poisson_draw = rpois(n(), ps_selection_prob),
        selected_ps = poisson_draw > 0
      )
    
    ps_set <- which(pop_data$selected_ps)
    cat("\nPS size from Poisson sampling:", length(ps_set), "\n")
    cat("Selection probability range:", 
        round(range(selection_prob), 3), "\n")
    
    assign_ps_weights <- function(data) {
      data %>% mutate(survey_weight = 1 / ps_selection_prob)
    }
  }
  
  # =====================================================
  # NON-PROBABILITY SAMPLE (NPS) SELECTION
  # =====================================================
  
  cat("\n--- NPS Selection ---\n")
  cat("NPS selection parameters:\n")
  cat("nps_betas:", round(nps_betas, 3), "\n")
  cat("nps_sigma:", nps_sigma, "\n")
  
  # Create design matrix WITH intercept (matching age_only function)
  X_nps_selection <- model.matrix(~ age_category, data = pop_data)
  
  # Validate dimensions
  if (ncol(X_nps_selection) != length(nps_betas)) {
    stop("Mismatch: design matrix has ", ncol(X_nps_selection), 
         " columns but nps_betas has length ", length(nps_betas))
  }
  
  # Generate errors and calculate selection probabilities
  nps_errors <- rnorm(nrow(pop_data), 0, nps_sigma)
  logit_p_nps <- X_nps_selection %*% nps_betas + nps_errors
  
  # Add MNAR effect if requested
  if (use_mnar) {
    cat("\nApplying MNAR mechanism based on Y_2\n")
    logit_p_nps <- logit_p_nps + mnar_shift * pop_data$Y_2
    cat("MNAR shift parameter:", mnar_shift, "\n")
  }
  
  nps_selection_prob <- plogis(logit_p_nps)
  
  # Report selection probabilities by age
  cat("\nNPS selection probabilities by age:\n")
  selection_summary <- pop_data %>%
    mutate(selection_prob = as.numeric(nps_selection_prob)) %>%
    group_by(age_category) %>%
    summarise(
      mean_prob = mean(selection_prob),
      min_prob = min(selection_prob),
      max_prob = max(selection_prob),
      .groups = "drop"
    )
  print(selection_summary)
  
  # If MNAR, also report by Y_2
  if (use_mnar) {
    cat("\nNPS selection probabilities by Y_2:\n")
    y2_summary <- pop_data %>%
      mutate(selection_prob = as.numeric(nps_selection_prob)) %>%
      group_by(Y_2) %>%
      summarise(
        mean_prob = mean(selection_prob),
        min_prob = min(selection_prob),
        max_prob = max(selection_prob),
        .groups = "drop"
      )
    print(y2_summary)
  }
  
  # =====================================================
  # MANAGE OVERLAP BETWEEN SAMPLES
  # =====================================================
  
  overlap_tolerance <- 0.02
  min_overlap_size <- round(length(ps_set) * (overlap_pct - overlap_tolerance))
  max_overlap_size <- round(length(ps_set) * (overlap_pct + overlap_tolerance))
  target_overlap_size <- round(length(ps_set) * overlap_pct)
  
  cat("\nTarget overlap:", target_overlap_size, 
      "acceptable range: [", min_overlap_size, ",", max_overlap_size, "]\n")
  
  # Sample NPS with overlap control
  max_iterations <- 100
  iteration <- 0
  overlap_found <- FALSE
  
  while (!overlap_found && iteration < max_iterations) {
    iteration <- iteration + 1
    
    nps_set <- sample(1:nrow(pop_data), N_nps, 
                      prob = nps_selection_prob, replace = FALSE)
    
    overlap <- intersect(ps_set, nps_set)
    overlap_size <- length(overlap)
    
    if (overlap_size >= min_overlap_size && overlap_size <= max_overlap_size) {
      overlap_found <- TRUE
      cat("Found acceptable overlap in iteration", iteration, 
          "- size:", overlap_size, "\n")
    }
  }
  
  if (!overlap_found) {
    cat("Warning: Could not achieve target overlap after", max_iterations, 
        "iterations. Using last sample.\n")
  }
  
  # =====================================================
  # CREATE FINAL DATASETS
  # =====================================================
  
  ps_only_set <- setdiff(ps_set, overlap)
  nps_only_set <- setdiff(nps_set, overlap)
  
  # Create datasets with appropriate missingness
  ps_data <- pop_data[ps_only_set, ] %>%
    mutate(Y_1 = NA)  # Y_1 missing in PS
  
  overlap_data <- pop_data[overlap, ]  # Both Y_1 and Y_2 observed
  
  nps_data <- pop_data[nps_only_set, ] %>%
    mutate(Y_2 = NA)  # Y_2 missing in NPS
  
  # Assign weights
  ps_data <- assign_ps_weights(ps_data)
  overlap_data <- assign_ps_weights(overlap_data)
  
  # For Poisson sampling, normalize weights
  if (sampling_scheme == "PO") {
    total_pop_size <- nrow(pop_data)
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
    current_sum <- sum(ps_with_overlap$survey_weight)
    
    cat("\nBefore normalization - Sum of survey weights:", current_sum, "\n")
    cat("Population size:", total_pop_size, "\n")
    
    normalization_factor <- total_pop_size / current_sum
    ps_data <- ps_data %>% mutate(survey_weight = survey_weight * normalization_factor)
    overlap_data <- overlap_data %>% mutate(survey_weight = survey_weight * normalization_factor)
    
    cat("Weight normalization factor:", round(normalization_factor, 3), "\n")
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
    cat("After normalization - Sum of survey weights:", sum(ps_with_overlap$survey_weight), "\n")
  } else {
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
  }
  
  # =====================================================
  # CREATE CELLS FROM PS+OVERLAP ONLY (WHERE Y_2 IS OBSERVED)
  # Key difference from age_only: cells include Y_2
  # =====================================================
  
  # Define cells based on (ST, age_category, Y_2)
  cell_summary <- ps_with_overlap %>%
    mutate(cell_id = interaction(ST, age_category, Y_2, drop = TRUE)) %>%
    group_by(cell_id) %>%
    summarise(
      n_ps_c = n(),
      w_bar_c = mean(survey_weight),
      ST = first(ST),
      age_category = first(age_category),
      Y_2 = first(Y_2),
      .groups = "drop"
    ) %>%
    mutate(cell_index = row_number())
  
  cat("\nNumber of cells from PS+overlap (ST × age × Y_2):", nrow(cell_summary), "\n")
  
  # Natural sort states
  natural_sort <- function(x) {
    if (all(grepl("^ST\\d+$", x))) {
      nums <- as.numeric(gsub("ST", "", x))
      x[order(nums)]
    } else {
      sort(x)  # For PUMA data
    }
  }
  
  unique_states <- natural_sort(unique(pop_data$ST))
  state_to_index <- setNames(1:length(unique_states), unique_states)
  
  # Convert to indices
  group_ps <- state_to_index[ps_data$ST]
  group_nps <- state_to_index[nps_data$ST]
  group_overlap <- state_to_index[overlap_data$ST]
  cell_to_group <- state_to_index[cell_summary$ST]
  
  # Create design matrices (age only - Y_2 is not included as covariate)
  X_ps <- model.matrix(~ age_category, data = ps_data)
  X_nps <- model.matrix(~ age_category, data = nps_data)
  X_overlap <- model.matrix(~ age_category, data = overlap_data)
  
  # Create cell design matrix
  X_cell <- model.matrix(~ age_category, data = cell_summary)
  
  # =====================================================
  # CREATE MINIMAL STAN DATA
  # =====================================================
  
  stan_data <- list(
    # Sample sizes
    N_ps = nrow(ps_data),
    N_nps = nrow(nps_data),
    N_overlap = nrow(overlap_data),
    J = length(unique_states),
    K = ncol(X_ps),
    
    # Group indicators
    group_ps = group_ps,
    group_nps = group_nps,
    group_overlap = group_overlap,
    
    # Design matrices
    X_ps = X_ps,
    X_nps = X_nps,
    X_overlap = X_overlap,
    
    # Observed outcomes
    Y2_ps = ps_data$Y_2,
    Y1_nps = nps_data$Y_1,
    Y1_overlap = overlap_data$Y_1,
    Y2_overlap = overlap_data$Y_2,
    
    # Cell information
    C = nrow(cell_summary),
    n_ps_c = cell_summary$n_ps_c,
    w_bar_c = cell_summary$w_bar_c,
    cell_to_group = cell_to_group,
    X_cell = X_cell
  )
  
  # =====================================================
  # SUMMARY REPORT
  # =====================================================
  
  cat("\n=== SAMPLING SUMMARY ===\n")
  cat("PS sample size:", nrow(ps_data), "\n")
  cat("NPS sample size:", nrow(nps_data), "\n")
  cat("Overlap size:", nrow(overlap_data), "\n")
  cat("Number of cells from PS+overlap:", nrow(cell_summary), "\n")
  cat("Number of states:", length(unique_states), "\n")
  cat("Number of covariates (K):", ncol(X_ps), "\n")
  cat("Design matrix colnames:", colnames(X_ps), "\n")
  
  # Analyze MNAR effect if applicable
  if (use_mnar) {
    cat("\n=== MNAR EFFECT ANALYSIS ===\n")
    y2_prop_pop <- mean(pop_data$Y_2)
    y2_prop_ps <- mean(ps_data$Y_2, na.rm = TRUE)
    y1_prop_nps <- mean(nps_data$Y_1, na.rm = TRUE)
    y2_prop_overlap <- mean(overlap_data$Y_2, na.rm = TRUE)
    
    cat("Y_2 proportion in population:", round(y2_prop_pop, 4), "\n")
    cat("Y_2 proportion in PS sample:", round(y2_prop_ps, 4), 
        "(diff:", round(y2_prop_ps - y2_prop_pop, 4), ")\n")
    cat("Y_1 proportion in NPS sample:", round(y1_prop_nps, 4), "\n")
    cat("Y_2 proportion in overlap:", round(y2_prop_overlap, 4), 
        "(diff:", round(y2_prop_overlap - y2_prop_pop, 4), ")\n")
  }
  
  cat("\nNote: NPS cell assignment will happen after Y_2 imputation in R\n")
  cat("      EMRP uses cells based on (ST × age × Y_2) for finer stratification\n")
  
  # Return results
  return(list(
    # Datasets
    ps_data = ps_data,
    nps_data = nps_data,
    overlap_data = overlap_data,
    
    # Cell information
    cell_summary = cell_summary,
    
    # Stan data
    stan_data = stan_data,
    
    # Combined PS data
    ps_with_overlap = ps_with_overlap,
    
    # Metadata
    actual_overlap_pct = length(overlap) / length(ps_set),
    unique_states = unique_states,
    state_to_index = state_to_index,
    sampling_scheme = sampling_scheme,
    use_mnar = use_mnar,
    
    # Selection parameters for reference
    ps_selection_params = list(
      scheme = sampling_scheme,
      betas = if(sampling_scheme == "PO") ps_betas else NULL,
      sigma = if(sampling_scheme == "PO") ps_sigma else NULL
    ),
    nps_selection_params = list(
      betas = nps_betas,
      sigma = nps_sigma,
      mnar_shift = if (use_mnar) mnar_shift else NULL
    )
  ))
}

samples_EMRP_predictive_ <- create_overlap_samples_EMRP_predictive(pop_data = pop_data_predictive, 
                                                                   sampling_scheme = "Stratified", use_mnar = F, seed = 17,
                                                                   overlap_pct = 0.11)

emrp_y2_predictor_model <- cmdstan_model("EMRP_y2_predictor.stan")


run_simulation_EMRP_predictive <- function(s, pop_data, overlap_pct, true_values, sampling_scheme,
                                           use_mnar = FALSE) {
  # -----------------------------------------------------------
  # 1. Create Samples and Prepare Stan Data
  # -----------------------------------------------------------
  samples <- create_overlap_samples_EMRP_predictive(
    pop_data = pop_data_predictive,
    N_ps = 1100,
    N_nps = 4300,
    overlap_pct = overlap_pct,
    sampling_scheme = sampling_scheme,
    use_mnar = use_mnar,
    seed = s + 17
  )
  
  stan_data <- samples$stan_data
  
  # -----------------------------------------------------------
  # 2. Fit the EMRP Predictive Stan Model
  # -----------------------------------------------------------
  fit <- emrp_y2_predictor_model$sample(
    data            = stan_data,
    chains          = 2,
    parallel_chains = 2,
    iter_warmup     = 500,
    iter_sampling   = 500,
    seed            = s + 42,
    refresh         = 0
  )
  
  # -----------------------------------------------------------
  # 3. Extract Posterior Draws
  # -----------------------------------------------------------
  # Extract Y2 imputations for NPS
  Y2_nps_draws <- fit$draws("Y2_nps_imputed", format = "df") %>% as_tibble()
  # Extract population cell counts
  N_c_draws <- fit$draws("N_c", format = "df") %>% as_tibble()
  # Extract model parameters for diagnostics
  gamma_draws <- fit$draws("gamma", format = "df") %>% as_tibble()
  beta_draws <- fit$draws("beta", format = "df") %>% as_tibble()
  alpha_draws <- fit$draws("alpha", format = "draws_matrix") %>% as.numeric()
  
  num_draws <- nrow(Y2_nps_draws)
  
  # -----------------------------------------------------------
  # 4. Process Each MCMC Draw to Create Blended Estimates
  # -----------------------------------------------------------
  draw_results <- lapply(seq_len(num_draws), function(m) {
    # Get Y2 imputations for this draw
    y2_imputed <- Y2_nps_draws[m, paste0("Y2_nps_imputed[", 1:nrow(samples$nps_data), "]")] %>% 
      unlist() %>% 
      unname()
    
    # Create cells for NPS with imputed Y2
    nps_with_y2 <- samples$nps_data %>%
      mutate(
        Y_2_imputed = y2_imputed,
        cell_id = interaction(ST, age_category, Y_2_imputed, drop = TRUE)
      )
    
    # Also add overlap data with its observed Y2
    overlap_with_cells <- samples$overlap_data %>%
      mutate(
        cell_id = interaction(ST, age_category, Y_2, drop = TRUE)
      )
    
    # Combine NPS and overlap for cell statistics
    combined_nps <- bind_rows(
      nps_with_y2 %>% select(cell_id, Y_1, ST, age_category),
      overlap_with_cells %>% select(cell_id, Y_1, ST, age_category)
    )
    
    # Calculate NPS cell statistics
    nps_cell_stats <- combined_nps %>%
      group_by(cell_id) %>%
      summarise(
        n_c_nps = n(),
        y_bar_c_nps = mean(Y_1),
        .groups = "drop"
      )
    
    # Get current draw's parameters for model predictions
    gamma_m <- gamma_draws[m, ] %>% select(starts_with("gamma")) %>% unlist()
    beta_m <- beta_draws[m, ] %>% select(starts_with("beta")) %>% unlist()
    alpha_m <- alpha_draws[m]
    
    # Also get state random effects for this draw
    u_draws <- fit$draws("u", format = "df") %>% as_tibble()
    v_draws <- fit$draws("v", format = "df") %>% as_tibble()
    u_m <- u_draws[m, ] %>% select(starts_with("u")) %>% unlist()
    v_m <- v_draws[m, ] %>% select(starts_with("v")) %>% unlist()
    
    # Map to cell structure and calculate blended estimates
    cell_data <- samples$cell_summary %>%
      mutate(cell_idx = row_number()) %>%
      left_join(nps_cell_stats, by = "cell_id") %>%
      mutate(
        n_c_nps = replace_na(n_c_nps, 0),
        y_bar_c_nps = replace_na(y_bar_c_nps, 0),
        N_c = N_c_draws[m, paste0("N_c[", cell_idx, "]")] %>% unlist() %>% unname()
      )
    
    # Calculate model predictions for each cell
    # Need to get state index for random effects
    state_indices <- samples$state_to_index[cell_data$ST]
    
    # Create design matrix for cells
    X_cell_pred <- model.matrix(~ age_category, data = cell_data)
    
    # Calculate predicted Y2 probability for each cell
    eta_Y2_cell <- X_cell_pred %*% gamma_m + v_m[state_indices]
    p_Y2_cell <- plogis(eta_Y2_cell)
    
    # For cells with Y_2 = 0
    eta_Y1_Y2_0 <- X_cell_pred %*% beta_m + alpha_m * 0 + u_m[state_indices]
    p_Y1_Y2_0 <- plogis(eta_Y1_Y2_0)
    
    # For cells with Y_2 = 1
    eta_Y1_Y2_1 <- X_cell_pred %*% beta_m + alpha_m * 1 + u_m[state_indices]
    p_Y1_Y2_1 <- plogis(eta_Y1_Y2_1)
    
    # Model prediction based on actual Y_2 value in the cell
    cell_data <- cell_data %>%
      mutate(
        model_pred = if_else(Y_2 == 1, p_Y1_Y2_1, p_Y1_Y2_0),
        # Blended estimate
        cell_estimate = if_else(
          n_c_nps > 0,
          n_c_nps * y_bar_c_nps + (N_c - n_c_nps) * model_pred,
          N_c * model_pred
        )
      )
    
    # Calculate group-level estimates
    group_estimates <- cell_data %>%
      group_by(ST) %>%
      summarise(
        numerator = sum(cell_estimate),
        denominator = sum(N_c),
        estimate = numerator / denominator,
        .groups = "drop"
      )
    
    # Calculate overall estimate
    overall_estimate <- sum(cell_data$cell_estimate) / sum(cell_data$N_c)
    
    return(list(
      group_estimates = group_estimates,
      overall_estimate = overall_estimate,
      alpha = alpha_m  # Store alpha for diagnostics
    ))
  })
  
  # -----------------------------------------------------------
  # 5. Aggregate Results Across Draws
  # -----------------------------------------------------------
  group_results <- list()
  
  # Process each state
  all_states <- sort(unique(samples$cell_summary$ST))
  
  for (st in all_states) {
    true_value <- true_values[[as.character(st)]]
    
    # Extract EMRP estimates for this state across draws
    estimates <- sapply(draw_results, function(x) {
      if (st %in% x$group_estimates$ST) {
        x$group_estimates$estimate[x$group_estimates$ST == st]
      } else {
        NA
      }
    })
    
    # Remove NA values
    estimates <- estimates[!is.na(estimates)]
    
    if (length(estimates) > 0) {
      # Calculate mean and variance
      hat_emrp <- mean(estimates)
      T_emrp <- var(estimates)
      
      # Calculate bias and MSE
      bias_emrp <- hat_emrp - true_value
      mse_emrp <- T_emrp + bias_emrp^2
    } else {
      hat_emrp <- T_emrp <- bias_emrp <- mse_emrp <- NA
    }
    
    # Store state-level results
    group_results[[as.character(st)]] <- list(
      hat_theta_emrp = hat_emrp,
      T_s_emrp = T_emrp,
      bias_emrp = bias_emrp,
      mse_emrp = mse_emrp
    )
  }
  
  # Process overall estimates
  true_value_overall <- true_values[["overall"]]
  overall_estimates <- sapply(draw_results, function(x) x$overall_estimate)
  
  # Calculate overall statistics
  hat_emrp_overall <- mean(overall_estimates)
  T_emrp_overall <- var(overall_estimates)
  bias_emrp_overall <- hat_emrp_overall - true_value_overall
  mse_emrp_overall <- T_emrp_overall + bias_emrp_overall^2
  
  # Store overall results
  group_results[["overall"]] <- list(
    hat_theta_emrp = hat_emrp_overall,
    T_s_emrp = T_emrp_overall,
    bias_emrp = bias_emrp_overall,
    mse_emrp = mse_emrp_overall
  )
  
  # Extract alpha values for diagnostics
  alpha_values <- sapply(draw_results, function(x) x$alpha)
  cat("Mean alpha (Y2->Y1 effect):", mean(alpha_values), "\n")
  
  return(group_results)
}

# Function to run a full batch of EMRP predictive simulations
run_EMRP_predictive_simulations <- function(S, pop_data, true_values, overlap_pct = 0.11, 
                                            sampling_scheme, use_mnar = FALSE) {
  # Compile the Stan model once
  cat("Compiling EMRP predictive Stan model...\n")
  emrp_y2_predictor_model <<- cmdstan_model("EMRP_y2_predictor.stan")
  
  # Set up parallel processing
  plan(multisession, workers = 3)  # Adjust based on available cores
  
  # Use progressr to track progress
  progressr::with_progress({
    p <- progressr::progressor(along = 1:S)
    
    # Run simulations in parallel
    simulation_results <- future_lapply(1:S, function(x) {
      p(sprintf("Simulation %d completed", x))  # Update progress bar
      run_simulation_EMRP_predictive(
        s = x,
        pop_data = pop_data,
        overlap_pct = overlap_pct,
        true_values = true_values,
        sampling_scheme = sampling_scheme,
        use_mnar = use_mnar
      )
    }, future.seed = TRUE)
  })
  
  # Clean up parallel workers
  plan(sequential)
  gc()
  
  # Determine MAR/MNAR suffix for file names
  mnar_suffix <- if(use_mnar) "MNAR" else "MAR"
  
  # Save results with appropriate naming
  results_filename <- paste0(sampling_scheme, "-EMRP-Predictive-", mnar_suffix, 
                             "-Results-", overlap_pct*100, "pct-", S, "s.rds")
  saveRDS(simulation_results, results_filename)
  cat("Results saved to:", results_filename, "\n")
  
  # Calculate metrics
  EMRP_metrics <- calculate_metrics_unified(
    simulation_results = simulation_results,
    true_values = true_values,
    method = "emrp"
  )
  
  # Save metrics with appropriate naming
  metrics_filename <- paste0(sampling_scheme, "-EMRP-Predictive-", mnar_suffix, 
                             "-Metrics-", overlap_pct*100, "pct-", S, "s.rds")
  saveRDS(EMRP_metrics, metrics_filename)
  cat("Metrics saved to:", metrics_filename, "\n")
  
  return(list(results = simulation_results, metrics = EMRP_metrics))
}

# Function to run all EMRP predictive simulations
run_all_EMRP_predictive_simulations <- function(S = 100, pop_data, true_values, overlap_pct = 0.11) {
  
  # Define simulation parameters
  sampling_schemes <- c("SRS", "Stratified", "PO")
  mnar_settings <- c(TRUE, FALSE)  # TRUE for MNAR, FALSE for MAR
  
  # Store all results
  all_results <- list()
  
  # Run all combinations
  for (scheme in sampling_schemes) {
    for (use_mnar in mnar_settings) {
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat("Running EMRP Predictive:", scheme, if(use_mnar) "MNAR" else "MAR", "\n")
      cat(paste(rep("=", 60), collapse = ""), "\n\n")
      
      result <- run_EMRP_predictive_simulations(
        S = S, 
        pop_data = pop_data, 
        true_values = true_values,
        use_mnar = use_mnar,
        overlap_pct = overlap_pct, 
        sampling_scheme = scheme
      )
      
      # Store result with descriptive key
      key <- paste0(scheme, "_", if(use_mnar) "MNAR" else "MAR")
      all_results[[key]] <- result
    }
  }
  
  return(all_results)
}

# Run all EMRP predictive simulations
all_results <- run_all_EMRP_predictive_simulations(
  S = 100, 
  pop_data = pop_data_predictive_age, 
  true_values = true_values_predictive_age, 
  overlap_pct = 0.11
)
