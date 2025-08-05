# ============================================
# UNIFIED DGP FOR MODELS 1-4
# ============================================
# This function generates data for all four models to pin-point MRP failure
# Perfectly matches enhanced_y2_dgp.R when model = 4
# Model 1: Y1 = X*beta + alpha*Y2
# Model 2: Y1 = X*beta + alpha*Y2 + gamma_j (marginal area effects)
# Model 3: Y1 = X*beta + (alpha + X*delta)*Y2 + gamma_j (marginal area effects)
# Model 4: Y1 = X*beta + (alpha + gamma_j + X*delta)*Y2 (area × Y2 interactions)

generate_unified_dgp <- function(
    mi_data, 
    model = 4,  # Which model (1-4)
    seed = 123,
    # Y2 effect parameters
    alpha = 1.5,              # Global Y2 effect (all models)
    sigma_gamma = 2.5,        # SD of area effects (models 2-4)
    # Interaction parameters
    demo_y2_interaction_sd = 1.0  # SD for demographic × Y2 (models 3 & 4)
) {
  
  set.seed(seed)
  
  # Validate inputs
  if (!model %in% 1:4) {
    stop("Model must be 1, 2, 3, or 4")
  }
  
  # Set parameters based on model
  if (model == 1) {
    use_area_effects <- FALSE
    use_demo_interactions <- FALSE
    model_desc <- "Model 1: Y1 = X*beta + alpha*Y2"
  } else if (model == 2) {
    use_area_effects <- TRUE
    use_demo_interactions <- FALSE
    model_desc <- "Model 2: Y1 = X*beta + alpha*Y2 + gamma_j"
  } else if (model == 3) {
    use_area_effects <- TRUE
    use_demo_interactions <- TRUE
    model_desc <- "Model 3: Y1 = X*beta + (alpha + X*delta)*Y2 + gamma_j"
  } else if (model == 4) {
    use_area_effects <- TRUE
    use_demo_interactions <- TRUE
    model_desc <- "Model 4: Y1 = X*beta + (alpha + gamma_j + X*delta)*Y2"
  }
  
  cat("\n=== GENERATING UNIFIED DGP ===\n")
  cat("Model:", model_desc, "\n")
  cat("Key features:\n")
  cat("- Area effects (gamma_j):", ifelse(use_area_effects, paste0("N(0, ", sigma_gamma^2, ")"), "None"), "\n")
  cat("- How area effects enter:", 
      ifelse(model == 1, "N/A", 
             ifelse(model %in% c(2, 3), "Marginal (additive)", 
                    "Interacting with Y2")), "\n")
  cat("- Demographic × Y2 interaction:", ifelse(use_demo_interactions, "Yes", "No"), "\n\n")
  
  # Helper functions
  inv_logit <- function(x) plogis(x)
  
  # Setup
  n_population <- nrow(mi_data)
  unique_states <- sort(unique(mi_data$ST))
  n_states <- length(unique_states)
  
  cat("Population size:", n_population, "\n")
  cat("Number of areas:", n_states, "\n\n")
  
  # =====================================================
  # CREATE DESIGN MATRICES
  # =====================================================
  
  # Create consistent design matrix for demographics
  X <- model.matrix(~ race_category + gender_binary + age_category, data = mi_data)
  K <- ncol(X)
  
  # Create area matrix for reference
  area_matrix <- model.matrix(~ ST - 1, data = mi_data)
  
  cat("\n=== DESIGN MATRICES ===\n")
  cat("Demographic matrix X:", dim(X), "\n")
  cat("Columns:", paste(colnames(X), collapse = ", "), "\n")
  cat("\nArea matrix:", dim(area_matrix), "\n")
  cat("Number of areas:", ncol(area_matrix), "\n")
  
  # =====================================================
  # PARAMETERS FOR Y2 MODEL
  # =====================================================
  
  cat("\n=== Y2 MODEL PARAMETERS ===\n")
  
  # Y2 model: logit(P(Y2=1)) = X*gamma
  gamma <- rnorm(K, 0, 0.5)
  gamma[1] <- -1.0  # Intercept
  names(gamma) <- colnames(X)
  
  cat("Y2 coefficients (gamma):\n")
  print(round(gamma, 3))
  
  # =====================================================
  # PARAMETERS FOR Y1 MODEL  
  # =====================================================
  
  cat("\n=== Y1 MODEL PARAMETERS ===\n")
  
  # Main effects of demographics on Y1
  beta <- rnorm(K, 0, 0.5)
  beta[1] <- -2.0  # Intercept
  names(beta) <- colnames(X)
  
  cat("\nMain demographic effects (beta):\n")
  print(round(beta, 3))
  
  # Global Y2 effect
  cat("\nGlobal Y2 effect (alpha):", alpha, "\n")
  
  # Demographic × Y2 interactions
  if (use_demo_interactions) {
    delta <- rnorm(K, 0, demo_y2_interaction_sd)
    delta[1] <- 0  # No interaction with intercept
    names(delta) <- paste0(colnames(X), ":Y2")
    
    cat("\nDemographic × Y2 interactions (delta):\n")
    print(round(delta, 3))
  } else {
    delta <- rep(0, K)
    names(delta) <- paste0(colnames(X), ":Y2")
    cat("\nNo demographic × Y2 interactions (delta = 0)\n")
  }
  
  # Area-specific Y2 effects (gamma_j)
  if (use_area_effects) {
    gamma_j <- rnorm(n_states, 0, sigma_gamma)
    names(gamma_j) <- unique_states
    
    cat("\nArea × Y2 effects (gamma_j) ~ N(0,", sigma_gamma^2, ")\n")
    cat("Summary of area-specific Y2 effects:\n")
    print(summary(gamma_j))
    
    # Show total Y2 effect range by area (alpha + gamma_j)
    total_y2_effect_by_area <- alpha + gamma_j
    cat("\nTotal Y2 effect by area (alpha + gamma_j):\n")
    cat("Range: [", round(min(total_y2_effect_by_area), 2), ",", 
        round(max(total_y2_effect_by_area), 2), "]\n")
  } else {
    gamma_j <- rep(0, n_states)
    names(gamma_j) <- unique_states
    cat("\nNo area effects (gamma_j = 0)\n")
  }
  
  # =====================================================
  # GENERATE Y2
  # =====================================================
  
  cat("\n=== GENERATING Y2 ===\n")
  
  # Linear predictor for Y2 (no area effects)
  eta_Y2 <- as.vector(X %*% gamma)
  
  # Generate Y2
  p_Y2 <- inv_logit(eta_Y2)
  Y_2 <- rbinom(n_population, 1, p_Y2)
  
  cat("Y2 prevalence:", round(mean(Y_2), 3), "\n")
  cat("Y2 prevalence by area range: [", 
      round(min(tapply(Y_2, mi_data$ST, mean)), 3), ",",
      round(max(tapply(Y_2, mi_data$ST, mean)), 3), "]\n")
  
  # =====================================================
  # GENERATE Y1 WITH INTERACTIONS
  # =====================================================
  
  cat("\n=== GENERATING Y1 WITH INTERACTIONS ===\n")
  
  # Calculate Y2 effects for each individual
  # Total Y2 effect = alpha + X*delta + area*gamma_j
  
  # 1. Main Y2 effect (alpha)
  y2_main_effect <- alpha * Y_2
  
  # 2. Demographic × Y2 interactions (X*delta*Y2)
  # Note: delta[1] = 0 (no intercept interaction)
  y2_demo_interaction <- as.vector((X %*% delta) * Y_2)
  
  # 3. Area × Y2 interactions
  # Map area-specific effects to individuals
  area_indices <- match(mi_data$ST, unique_states)
  if (model == 4) {
    # Model 4: area effects interact with Y2
    y2_area_interaction <- gamma_j[area_indices] * Y_2
  } else {
    # Models 2-3: area effects are marginal (not interacting with Y2)
    y2_area_interaction <- rep(0, n_population)
  }
  
  # Total Y2 contribution
  total_y2_contribution <- y2_main_effect + y2_demo_interaction + y2_area_interaction
  
  # Linear predictor for Y1
  eta_Y1 <- as.vector(X %*% beta) + total_y2_contribution
  
  # Add marginal area effects for models 2 and 3
  if (model %in% c(2, 3)) {
    eta_Y1 <- eta_Y1 + gamma_j[area_indices]
  }
  
  # Generate Y1
  p_Y1 <- inv_logit(eta_Y1)
  Y_1 <- rbinom(n_population, 1, p_Y1)
  
  cat("Y1 prevalence:", round(mean(Y_1), 3), "\n")
  cat("Y1 prevalence by area range: [", 
      round(min(tapply(Y_1, mi_data$ST, mean)), 3), ",",
      round(max(tapply(Y_1, mi_data$ST, mean)), 3), "]\n")
  
  # =====================================================
  # CREATE OUTPUT DATA
  # =====================================================
  
  mi_data_with_outcomes <- mi_data %>%
    mutate(Y_1 = Y_1, Y_2 = Y_2)
  
  # Calculate true values by area
  true_values <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(mean_Y1 = mean(Y_1), .groups = "drop") %>%
    deframe()
  
  true_values <- c(true_values, overall = mean(Y_1))
  
  # =====================================================
  # ANALYZE EXPECTED BIAS PATTERNS
  # =====================================================
  
  cat("\n=== EXPECTED BIAS ANALYSIS ===\n")
  
  # Area-level analysis
  area_analysis <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(
      n = n(),
      Y1_prev = mean(Y_1),
      Y2_prev = mean(Y_2),
      # Demographics
      prop_female = mean(gender_binary == "Female"),
      prop_age_60plus = mean(age_category %in% c("60-74", "75+")),
      .groups = "drop"
    )
  
  # Add area effects and expected bias patterns based on model
  if (model == 1) {
    area_analysis <- area_analysis %>%
      mutate(
        gamma_j = 0,
        total_y2_effect = alpha,
        expected_bias_main = alpha * Y2_prev,
        demo_interaction_bias = 0
      )
  } else if (model == 2) {
    area_analysis <- area_analysis %>%
      mutate(
        gamma_j = gamma_j[ST],
        total_y2_effect = alpha,  # Area effects don't interact with Y2
        expected_bias_main = alpha * Y2_prev,
        demo_interaction_bias = 0
      )
  } else if (model == 3) {
    # Need to calculate average demographic interaction effect by area
    demo_bias_by_area <- numeric(nrow(area_analysis))
    for (i in 1:nrow(area_analysis)) {
      area <- area_analysis$ST[i]
      area_data <- mi_data_with_outcomes %>% filter(ST == area)
      X_area <- model.matrix(~ race_category + gender_binary + age_category, data = area_data)
      avg_demo_effect <- mean(X_area %*% delta)
      demo_bias_by_area[i] <- avg_demo_effect * area_analysis$Y2_prev[i]
    }
    
    area_analysis <- area_analysis %>%
      mutate(
        gamma_j = gamma_j[ST],
        total_y2_effect = alpha,  # Area effects don't interact with Y2
        expected_bias_main = alpha * Y2_prev,
        demo_interaction_bias = demo_bias_by_area
      )
  } else if (model == 4) {
    # Calculate average total Y2 effect for each area
    total_y2_effects <- numeric(nrow(area_analysis))
    demo_bias_by_area <- numeric(nrow(area_analysis))
    
    for (i in 1:nrow(area_analysis)) {
      area <- area_analysis$ST[i]
      area_data <- mi_data_with_outcomes %>% filter(ST == area)
      X_area <- model.matrix(~ race_category + gender_binary + age_category, data = area_data)
      
      # Total Y2 effect for each individual
      individual_effects <- alpha + gamma_j[area] + as.vector(X_area %*% delta)
      total_y2_effects[i] <- mean(individual_effects)
      
      # Demo interaction bias
      demo_bias_by_area[i] <- mean(X_area %*% delta) * area_analysis$Y2_prev[i]
    }
    
    area_analysis <- area_analysis %>%
      left_join(
        tibble(
          ST = names(gamma_j),
          gamma_j = gamma_j,
          total_y2_effect = total_y2_effects
        ),
        by = "ST"
      ) %>%
      mutate(
        expected_bias_main = total_y2_effect * Y2_prev,
        demo_interaction_bias = demo_bias_by_area
      )
  }
  
  cat("\nMRP will miss:\n")
  if (model == 1) {
    cat("- The Y2 effect (alpha =", alpha, ")\n")
  } else if (model == 2) {
    cat("- The Y2 effect (alpha =", alpha, ")\n")
    cat("- But can model area effects as marginal effects\n")
  } else if (model == 3) {
    cat("- The Y2 effect (alpha =", alpha, ")\n")
    cat("- Demographic × Y2 interactions\n")
    cat("- But can model area effects as marginal effects\n")
  } else if (model == 4) {
    cat("- ALL Y2-related effects:\n")
    cat("  1. Main Y2 effect (alpha =", alpha, ")\n")
    cat("  2. Area × Y2 interactions (SD =", sigma_gamma, ")\n")
    cat("  3. Demographic × Y2 interactions\n")
    cat("- This will create systematic area-level bias in MRP estimates.\n")
  }
  
  # Show correlation between Y2 prevalence and Y1
  cor_y1_y2 <- cor(area_analysis$Y1_prev, area_analysis$Y2_prev)
  cat("\nArea-level correlation between Y1 and Y2 prevalence:", round(cor_y1_y2, 3), "\n")
  
  # =====================================================
  # STORE ALL PARAMETERS
  # =====================================================
  
  true_params <- list(
    # Model info
    model = model,
    model_desc = model_desc,
    
    # Y2 model parameters
    gamma = gamma,
    sigma_v = 0,  # 0 by design
    
    # Y1 model parameters
    beta = beta,
    alpha = alpha,
    delta = delta,  # Demographic × Y2 interactions
    gamma_j = gamma_j,  # Area × Y2 interactions
    sigma_gamma = sigma_gamma,
    sigma_u = ifelse(model %in% c(2, 3), sigma_gamma, 0),  # 0 by design
    
    # Total effects by area
    total_y2_effect = if (use_area_effects) total_y2_effect_by_area else rep(alpha, n_states),
    
    # For diagnostics
    area_analysis = area_analysis,
    
    # Design info
    K = K,
    J = n_states
  )
  
  # =====================================================
  # FINAL SUMMARY
  # =====================================================
  
  cat("\n=== SUMMARY ===\n")
  cat("Model", model, "DGP created:\n")
  cat("- Y2 prevalence:", round(mean(Y_2), 3), "\n")
  cat("- Y1 prevalence:", round(mean(Y_1), 3), "\n")
  cat("- Correlation(Y1, Y2):", round(cor(Y_1, Y_2), 3), "\n")
  
  if (model == 1) {
    cat("\nExpected outcome:\n")
    cat("- Both OG and MRP should perform well\n")
    cat("- This is the baseline case\n")
  } else if (model == 2) {
    cat("\nExpected outcome:\n")
    cat("- MRP can model area effects (gamma_j as marginal effects)\n")
    cat("- Both methods should perform reasonably well\n")
  } else if (model == 3) {
    cat("\nExpected outcome:\n")
    cat("- MRP can model area effects but misses demographic × Y2 interactions\n")
    cat("- OG should show some advantage\n")
  } else if (model == 4) {
    cat("\nExpected outcome:\n")
    cat("- MRP cannot capture area × Y2 interactions\n")
    cat("- OG should show substantial advantage\n")
    cat("- MRP coverage should fail dramatically\n")
  }
  
  return(list(
    data = mi_data_with_outcomes,
    true_values = true_values,
    parameters = true_params,
    design_matrices = list(X = X, area_matrix = area_matrix)
  ))
}

# ============================================
# EXAMPLE USAGE
# ============================================

if (FALSE) {  # Set to TRUE to run examples
  
  # Load Michigan data
  # mi_analysis_data <- readRDS("mi_pums_analysis_data.rds")
  
  # Generate all four models
  for (model_num in 1:4) {
    cat("\n", rep("=", 60), "\n")
    cat("GENERATING MODEL", model_num, "\n")
    cat(rep("=", 60), "\n")
    
    result <- generate_unified_dgp(
      mi_data = mi_analysis_data,
      model = model_num,
      # seed = 2024 + model_num,
      seed = 17,
      alpha = 1.5,
      sigma_gamma = 2.5,
      demo_y2_interaction_sd = 1.0
    )
    
    # Save results
    saveRDS(result$data, paste0("pop_data_model_", model_num, ".rds"))
    saveRDS(result$true_values, paste0("true_values_model_", model_num, ".rds"))
    saveRDS(result$parameters, paste0("true_params_model_", model_num, ".rds"))
  }
}
