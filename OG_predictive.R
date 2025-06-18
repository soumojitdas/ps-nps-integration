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

create_overlap_samples_age_only <- function(pop_data,
                                            N_ps = 1100,
                                            N_nps = 4300,
                                            overlap_pct = 0.11,
                                            sampling_scheme = "SRS",  # Options: "SRS", "Stratified", "PO"
                                            # PS selection parameters (full vector with intercept first)
                                            ps_betas = c(-3.1, 0.2, 0.3, 0.4, 0.2),  # length 5: intercept + 4 age effects
                                            ps_sigma = 0.1,
                                            # NPS selection parameters (full vector with intercept first)
                                            nps_betas = c(-2.703002, 0.1404789, 0.4330459, 0.4223657, 0.2616391),  # length 5
                                            nps_sigma = 0.1,
                                            # MNAR mechanism parameters
                                            use_mnar = FALSE,
                                            mnar_shift = 0.5,   # Simple shift in selection prob based on Y_2
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
  
  cat("=== AGE-ONLY SAMPLING (ST × age_category cells) ===\n")
  cat("Sampling scheme:", sampling_scheme, "\n")
  cat("MNAR mechanism:", ifelse(use_mnar, "ENABLED", "DISABLED"), "\n")
  
  # Create cells based on ST × age_category only (KEY DIFFERENCE FROM EMRP)
  poststrat_cells <- pop_data %>%
    mutate(cell_id = interaction(ST, age_category, drop = TRUE)) %>%
    group_by(cell_id, ST, age_category) %>%
    summarise(
      n_pop = n(),
      Y_1_mean = mean(Y_1),
      Y_2_mean = mean(Y_2),
      .groups = "drop"
    ) %>%
    mutate(cell_index = as.integer(factor(cell_id)))
  
  cat("\nNumber of cells (ST × age):", nrow(poststrat_cells), "\n")
  
  # Create cell-level design matrix WITH INTERCEPT
  X_cell <- model.matrix(~ age_category, data = poststrat_cells)
  
  # Map population to cells
  pop_data <- pop_data %>%
    mutate(cell_id = interaction(ST, age_category, drop = TRUE)) %>%
    left_join(poststrat_cells %>% select(cell_id, cell_index, n_pop), by = "cell_id")
  
  # =====================================================
  # PROBABILITY SAMPLE (PS) SELECTION
  # =====================================================
  
  if (sampling_scheme == "SRS") {
    # Simple Random Sampling
    ps_set <- sample(1:nrow(pop_data), N_ps, replace = FALSE)
    cat("\nPS size from Simple Random Sampling:", length(ps_set), "\n")
    
    # Constant weight for SRS
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
    
    # Generate individual-level errors
    ps_errors <- rnorm(nrow(pop_data), 0, ps_sigma)
    
    # Calculate selection probabilities
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
  
  # Create design matrix WITH intercept
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
  # MANAGE OVERLAP
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
  
  # For Poisson sampling, normalize weights to sum to population size
  if (sampling_scheme == "PO") {
    total_pop_size <- nrow(pop_data)
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
    current_sum <- sum(ps_with_overlap$survey_weight)
    
    cat("\nBefore normalization - Sum of survey weights:", current_sum, "\n")
    cat("Population size:", total_pop_size, "\n")
    
    # Apply normalization factor
    normalization_factor <- total_pop_size / current_sum
    ps_data <- ps_data %>% mutate(survey_weight = survey_weight * normalization_factor)
    overlap_data <- overlap_data %>% mutate(survey_weight = survey_weight * normalization_factor)
    
    cat("Weight normalization factor:", round(normalization_factor, 3), "\n")
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
    cat("After normalization - Sum of survey weights:", sum(ps_with_overlap$survey_weight), "\n")
  } else {
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
  }
  
  # Calculate cell-level statistics
  cell_summary <- ps_with_overlap %>%
    group_by(cell_index) %>%
    summarise(
      n_ps_c = n(),
      w_bar_c = mean(survey_weight),
      cell_id = first(cell_id),
      ST = first(ST),
      .groups = "drop"
    )
  
  cell_summary <- left_join(cell_summary, 
                            poststrat_cells |> dplyr::select(cell_index, n_pop),
                            by = c("cell_index"))
  
  # Natural sort function
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
  
  # Create design matrices - ALL WITH INTERCEPT
  X_ps <- model.matrix(~ age_category, data = ps_data)
  X_nps <- model.matrix(~ age_category, data = nps_data)
  X_overlap <- model.matrix(~ age_category, data = overlap_data)
  
  # Align X_cell matrix
  X_cell_aligned <- matrix(0, nrow = nrow(cell_summary), ncol = ncol(X_cell))
  colnames(X_cell_aligned) <- colnames(X_cell)
  
  for (i in 1:nrow(cell_summary)) {
    cell_idx <- cell_summary$cell_index[i]
    poststrat_row <- which(poststrat_cells$cell_index == cell_idx)
    if (length(poststrat_row) > 0) {
      X_cell_aligned[i, ] <- X_cell[poststrat_row, ]
    }
  }
  
  # Map cells for NPS and overlap
  nps_cell_id <- match(nps_data$cell_index, cell_summary$cell_index)
  overlap_cell_id <- match(overlap_data$cell_index, cell_summary$cell_index)
  
  # Filter observations without cell mapping
  valid_nps <- !is.na(nps_cell_id)
  valid_overlap <- !is.na(overlap_cell_id)
  
  if (any(!valid_nps)) {
    cat("\nFiltering", sum(!valid_nps), "NPS observations without PS cell mapping\n")
  }
  if (any(!valid_overlap)) {
    cat("Filtering", sum(!valid_overlap), "overlap observations without PS cell mapping\n")
  }
  
  # Create complete Stan data
  stan_data <- list(
    N_ps = nrow(ps_data),
    N_nps = sum(valid_nps),
    N_overlap = sum(valid_overlap),
    J = length(unique_states),
    K = ncol(X_ps),  # Should be 5 (intercept + 4 age dummies)
    group_ps = group_ps,
    group_nps = group_nps[valid_nps],
    group_overlap = group_overlap[valid_overlap],
    X_ps = X_ps,
    X_nps = X_nps[valid_nps, , drop = FALSE],
    X_overlap = X_overlap[valid_overlap, , drop = FALSE],
    Y2_ps = ps_data$Y_2,
    Y1_nps = nps_data$Y_1[valid_nps],
    Y1_overlap = overlap_data$Y_1[valid_overlap],
    Y2_overlap = overlap_data$Y_2[valid_overlap],
    C = nrow(cell_summary),
    N_c = cell_summary$n_pop,  # Known N_c for diagnostics
    n_ps_c = cell_summary$n_ps_c,
    w_bar_c = cell_summary$w_bar_c,
    cell_to_group = cell_to_group,
    X_cell = X_cell_aligned,
    nps_cell_id = nps_cell_id[valid_nps],
    overlap_cell_id = overlap_cell_id[valid_overlap],
    ps_survey_weights = ps_data$survey_weight,
    overlap_survey_weights = overlap_data$survey_weight[valid_overlap]
  )
  
  # =====================================================
  # SUMMARY REPORT
  # =====================================================
  
  cat("\n=== SAMPLING SUMMARY ===\n")
  cat("PS sample size:", nrow(ps_data), "\n")
  cat("NPS sample size:", nrow(nps_data), "(", sum(valid_nps), "after filtering)\n")
  cat("Overlap size:", nrow(overlap_data), "(", sum(valid_overlap), "after filtering)\n")
  cat("Actual overlap percentage:", round(100 * length(overlap) / length(ps_set), 2), "%\n")
  cat("Number of cells in PS:", nrow(cell_summary), "\n")
  cat("Number of covariates (K):", ncol(X_ps), "\n")
  cat("Design matrix colnames:", colnames(X_ps), "\n")
  cat("Number of small areas:", length(unique_states), "\n")
  
  # Analyze MNAR effect on sample composition if applicable
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
  
  # Return everything needed
  return(list(
    # Datasets
    ps_data = ps_data,
    nps_data = nps_data[valid_nps, ],
    overlap_data = overlap_data[valid_overlap, ],
    
    # Cell information
    poststrat_cells = poststrat_cells,
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

{
# possibly, the corrected one ----

create_overlap_samples_age_only <- function(pop_data,
                                            N_ps = 1100,
                                            N_nps = 4300,
                                            overlap_pct = 0.11,
                                            sampling_scheme = "SRS",  # Options: "SRS", "Stratified", "PO"
                                            # PS selection parameters (full vector with intercept first)
                                            ps_betas = c(-3.1, 0.2, 0.3, 0.4, 0.2),  # length 5: intercept + 4 age effects
                                            ps_sigma = 0.1,
                                            # NPS selection parameters (full vector with intercept first)
                                            nps_betas = c(-2.703002, 0.1404789, 0.4330459, 0.4223657, 0.2616391),  # length 5
                                            nps_sigma = 0.1,
                                            # MNAR mechanism parameters
                                            use_mnar = FALSE,
                                            mnar_shift = 0.5,   # Simple shift in selection prob based on Y_2
                                            seed = 17) {
  set.seed(seed)

  # # Validate inputs
  # required_vars <- c("Y_1", "Y_2", "ST", "age_category", "race_category", "gender_binary")
  # missing_vars <- setdiff(required_vars, names(pop_data))
  # if (length(missing_vars) > 0) {
  #   stop("Population data missing required variables: ", paste(missing_vars, collapse = ", "))
  # }
  #
  # # Validate beta vectors
  # if (length(ps_betas) != 5) {
  #   stop("ps_betas must have length 5 (intercept + 4 age effects)")
  # }
  # if (length(nps_betas) != 5) {
  #   stop("nps_betas must have length 5 (intercept + 4 age effects)")
  # }

  # Verify sampling scheme
  sampling_scheme <- match.arg(sampling_scheme, c("SRS", "Stratified", "PO"))

  # Ensure age factor is properly set
  pop_data <- pop_data %>%
    mutate(
      age_category = factor(age_category,
                            levels = c("18-29", "30-44", "45-59", "60-74", "75+"))
    )

  cat("=== AGE-ONLY SAMPLING with FULL MODEL ===\n")
  cat("Selection based on: age_category only\n")
  cat("Model/cells based on: ST × race × gender × age\n")
  cat("Sampling scheme:", sampling_scheme, "\n")
  cat("MNAR mechanism:", ifelse(use_mnar, "ENABLED", "DISABLED"), "\n")

  # Create cells based on ST × race_category × gender_binary × age_category (FULL STRATIFICATION)
  poststrat_cells <- pop_data %>%
    mutate(cell_id = interaction(ST, race_category, gender_binary, age_category, drop = TRUE)) %>%
    group_by(cell_id, ST, race_category, gender_binary, age_category) %>%
    summarise(
      n_pop = n(),
      Y_1_mean = mean(Y_1),
      Y_2_mean = mean(Y_2),
      .groups = "drop"
    ) %>%
    mutate(cell_index = as.integer(factor(cell_id)))

  cat("\nNumber of cells (ST × race × gender × age):", nrow(poststrat_cells), "\n")

  # Create cell-level design matrix WITH FULL MODEL
  X_cell <- model.matrix(~ race_category + gender_binary + age_category, data = poststrat_cells)

  # Map population to cells
  pop_data <- pop_data %>%
    mutate(cell_id = interaction(ST, race_category, gender_binary, age_category, drop = TRUE)) %>%
    left_join(poststrat_cells %>% select(cell_id, cell_index, n_pop), by = "cell_id")

  # =====================================================
  # PROBABILITY SAMPLE (PS) SELECTION
  # =====================================================

  if (sampling_scheme == "SRS") {
    # Simple Random Sampling
    ps_set <- sample(1:nrow(pop_data), N_ps, replace = FALSE)
    cat("\nPS size from Simple Random Sampling:", length(ps_set), "\n")

    # Constant weight for SRS
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
    # Poisson sampling using age ONLY for selection
    cat("\nPS Poisson sampling parameters:\n")
    cat("ps_betas:", round(ps_betas, 3), "\n")
    cat("ps_sigma:", ps_sigma, "\n")

    # Create design matrix for SELECTION (age only)
    X_ps_selection <- model.matrix(~ age_category, data = pop_data)

    # Validate dimensions
    if (ncol(X_ps_selection) != length(ps_betas)) {
      stop("Mismatch: design matrix has ", ncol(X_ps_selection),
           " columns but ps_betas has length ", length(ps_betas))
    }

    # Generate individual-level errors
    ps_errors <- rnorm(nrow(pop_data), 0, ps_sigma)

    # Calculate selection probabilities
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

  # Create design matrix for SELECTION (age only)
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
  # MANAGE OVERLAP
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

  # For Poisson sampling, normalize weights to sum to population size
  if (sampling_scheme == "PO") {
    total_pop_size <- nrow(pop_data)
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
    current_sum <- sum(ps_with_overlap$survey_weight)

    cat("\nBefore normalization - Sum of survey weights:", current_sum, "\n")
    cat("Population size:", total_pop_size, "\n")

    # Apply normalization factor
    normalization_factor <- total_pop_size / current_sum
    ps_data <- ps_data %>% mutate(survey_weight = survey_weight * normalization_factor)
    overlap_data <- overlap_data %>% mutate(survey_weight = survey_weight * normalization_factor)

    cat("Weight normalization factor:", round(normalization_factor, 3), "\n")
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
    cat("After normalization - Sum of survey weights:", sum(ps_with_overlap$survey_weight), "\n")
  } else {
    ps_with_overlap <- bind_rows(ps_data, overlap_data)
  }

  # Calculate cell-level statistics
  cell_summary <- ps_with_overlap %>%
    group_by(cell_index) %>%
    summarise(
      n_ps_c = n(),
      w_bar_c = mean(survey_weight),
      cell_id = first(cell_id),
      ST = first(ST),
      .groups = "drop"
    )

  cell_summary <- left_join(cell_summary,
                            poststrat_cells |> dplyr::select(cell_index, n_pop),
                            by = c("cell_index"))

  # Natural sort function
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

  # Create design matrices - NOW WITH FULL COVARIATES
  X_ps <- model.matrix(~ race_category + gender_binary + age_category, data = ps_data)
  X_nps <- model.matrix(~ race_category + gender_binary + age_category, data = nps_data)
  X_overlap <- model.matrix(~ race_category + gender_binary + age_category, data = overlap_data)

  # Align X_cell matrix
  X_cell_aligned <- matrix(0, nrow = nrow(cell_summary), ncol = ncol(X_cell))
  colnames(X_cell_aligned) <- colnames(X_cell)

  for (i in 1:nrow(cell_summary)) {
    cell_idx <- cell_summary$cell_index[i]
    poststrat_row <- which(poststrat_cells$cell_index == cell_idx)
    if (length(poststrat_row) > 0) {
      X_cell_aligned[i, ] <- X_cell[poststrat_row, ]
    }
  }

  # Map cells for NPS and overlap
  nps_cell_id <- match(nps_data$cell_index, cell_summary$cell_index)
  overlap_cell_id <- match(overlap_data$cell_index, cell_summary$cell_index)

  # Filter observations without cell mapping
  valid_nps <- !is.na(nps_cell_id)
  valid_overlap <- !is.na(overlap_cell_id)

  if (any(!valid_nps)) {
    cat("\nFiltering", sum(!valid_nps), "NPS observations without PS cell mapping\n")
  }
  if (any(!valid_overlap)) {
    cat("Filtering", sum(!valid_overlap), "overlap observations without PS cell mapping\n")
  }

  # Create complete Stan data
  stan_data <- list(
    N_ps = nrow(ps_data),
    N_nps = sum(valid_nps),
    N_overlap = sum(valid_overlap),
    J = length(unique_states),
    K = ncol(X_ps),  # Now includes race, gender, and age
    group_ps = group_ps,
    group_nps = group_nps[valid_nps],
    group_overlap = group_overlap[valid_overlap],
    X_ps = X_ps,
    X_nps = X_nps[valid_nps, , drop = FALSE],
    X_overlap = X_overlap[valid_overlap, , drop = FALSE],
    Y2_ps = ps_data$Y_2,
    Y1_nps = nps_data$Y_1[valid_nps],
    Y1_overlap = overlap_data$Y_1[valid_overlap],
    Y2_overlap = overlap_data$Y_2[valid_overlap],
    C = nrow(cell_summary),
    N_c = cell_summary$n_pop,  # Known N_c for diagnostics
    n_ps_c = cell_summary$n_ps_c,
    w_bar_c = cell_summary$w_bar_c,
    cell_to_group = cell_to_group,
    X_cell = X_cell_aligned,
    nps_cell_id = nps_cell_id[valid_nps],
    overlap_cell_id = overlap_cell_id[valid_overlap],
    ps_survey_weights = ps_data$survey_weight,
    overlap_survey_weights = overlap_data$survey_weight[valid_overlap]
  )

  # =====================================================
  # SUMMARY REPORT
  # =====================================================

  cat("\n=== SAMPLING SUMMARY ===\n")
  cat("PS sample size:", nrow(ps_data), "\n")
  cat("NPS sample size:", nrow(nps_data), "(", sum(valid_nps), "after filtering)\n")
  cat("Overlap size:", nrow(overlap_data), "(", sum(valid_overlap), "after filtering)\n")
  cat("Actual overlap percentage:", round(100 * length(overlap) / length(ps_set), 2), "%\n")
  cat("Number of cells in PS:", nrow(cell_summary), "\n")
  cat("Number of covariates (K):", ncol(X_ps), "\n")
  cat("Design matrix colnames:", colnames(X_ps), "\n")
  cat("Number of small areas:", length(unique_states), "\n")

  # Analyze MNAR effect on sample composition if applicable
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

  # Return everything needed
  return(list(
    # Datasets
    ps_data = ps_data,
    nps_data = nps_data[valid_nps, ],
    overlap_data = overlap_data[valid_overlap, ],

    # Cell information
    poststrat_cells = poststrat_cells,
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
}
# ----

samples_OG_predictive_ <- create_overlap_samples_age_only(pop_data = pop_data_predictive, 
                                                          sampling_scheme = "SRS", use_mnar = F, seed = 17,
                                                          overlap_pct = 0.11)

y2_covariate_model <- cmdstan_model("stan_model_y2_covariate.stan")

# Reiter's rule ----
run_simulation_theta1_modified <- function(s, pop_data, overlap_pct, true_values, sampling_scheme,
                                           use_mnar = FALSE, use_mse = FALSE) {
  # -----------------------------------------------------------
  # 1. Create Samples and Prepare Stan Data
  # -----------------------------------------------------------
  # 1. Create samples with the complete data structure for Stan
  samples <- create_overlap_samples_age_only(
    pop_data = pop_data,
    N_ps = 1100,
    N_nps = 4300,
    overlap_pct = overlap_pct,
    sampling_scheme = sampling_scheme,
    use_mnar = use_mnar,
    seed = s + 17
  )
  
  ps_sample      <- samples$ps_data           # PS-only sample
  overlap_sample <- samples$overlap_data      # Overlap sample
  stan_data      <- samples$stan_data
  
  # -----------------------------------------------------------
  # 2. Fit the Modified Stan Model that returns probabilities
  # -----------------------------------------------------------
  fit <- y2_covariate_model$sample(
    data            = stan_data,
    chains          = 2,
    parallel_chains = 2,
    iter_warmup     = 500,
    iter_sampling   = 500,
    seed            = s + 42,
    refresh         = 100
  )
  
  # fit$time()
  # fit$diagnostic_summary()
  # print(fit$summary(c("gamma", "sigma_v", "beta", "sigma_u", "alpha")), n = 51)
  # true_params_predictive
  
  # bayesplot::mcmc_recover_hist(fit$draws(c("gamma", "sigma_v", "beta", "sigma_u", "alpha")),
  #                              unlist(true_params_predictive[c("gamma", "sigma_v", "beta", "sigma_u", "alpha")]))
  
  # # various PPC checks ----
  # # Extract the true Y_1 values with grouping variables
  # y1_data <- left_join(
  #   ps_sample |> select(person_id, ST, age_category, race_category, gender_binary),
  #   pop_data_predictive |> select(person_id, Y_1),
  #   by = "person_id"
  # )
  # 
  # # Get just the Y_1 values for basic PPCs
  # y1_true <- y1_data$Y_1
  # 
  # # Extract posterior predictive draws
  # y1_rep <- fit$draws("Y1_ps_imputed", format = "matrix")
  # 
  # # Basic PPC
  # bayesplot::ppc_bars(y1_true, y1_rep)
  # 
  # # Grouped PPCs by different variables:
  # 
  # # 1. By age category
  # bayesplot::ppc_bars_grouped(
  #   y = y1_true,
  #   yrep = y1_rep,
  #   group = y1_data$age_category
  # )
  # 
  # # 2. By race category
  # bayesplot::ppc_bars_grouped(
  #   y = y1_true,
  #   yrep = y1_rep,
  #   group = y1_data$race_category
  # )
  # 
  # # 3. By gender
  # bayesplot::ppc_bars_grouped(
  #   y = y1_true,
  #   yrep = y1_rep,
  #   group = y1_data$gender_binary
  # )
  # 
  # # 4. By state (might be too many groups - select a subset if needed)
  # bayesplot::ppc_bars_grouped(
  #   y = y1_true,
  #   yrep = y1_rep,
  #   group = y1_data$ST
  # )
  # 
  # # Check proportion of Y_1=1 by group
  # bayesplot::ppc_stat_grouped(
  #   y = y1_true,
  #   yrep = y1_rep,
  #   group = y1_data$age_category,
  #   stat = function(y) mean(y)  # proportion of 1s
  # )
  # 
  # # Multiple groups - create combined factor
  # y1_data <- y1_data |>
  #   mutate(race_gender_age = interaction(race_category, gender_binary, age_category))
  # 
  # bayesplot::ppc_bars_grouped(
  #   y = y1_true,
  #   yrep = y1_rep,
  #   group = y1_data$race_gender_age
  # )
  # 
  # 
  # 
  # -----------------------------------------------------------
  # 3. Extract Posterior Probability Draws (instead of binary imputations)
  # -----------------------------------------------------------
  # Extract probabilities P(Y1=1|Y2) instead of binary imputations
  posterior_probs <- fit$draws("prob_Y1_given_Y2", format = "df") |> as_tibble()
  posterior_mrp <- fit$draws("group_mrp_estimate", format = "df") |> as_tibble()
  overall_mrp <- fit$draws("overall_mrp_estimate", format = "draws_matrix") |> as.numeric()
  
  all_states <- sort(unique(pop_data$ST))
  num_draws  <- nrow(posterior_probs)
  group_results <- list()
  
  # -----------------------------------------------------------
  # 4. Calculate PS Design-Based Estimates for Each Draw
  # -----------------------------------------------------------
  group_ps_estimates <- lapply(seq_len(num_draws), function(m) {
    # Get imputed PROBABILITIES for PS-only sample in this draw
    Y_probs <- posterior_probs[m, 1:nrow(ps_sample)] |> unlist()
    ps_with_probs <- ps_sample |> mutate(Y_val = Y_probs)
    
    # For overlap sample, use the actual observed Y1 BINARY values
    # This preserves the actual observed outcomes where available
    overlap_with_Y1 <- overlap_sample |> mutate(Y_val = as.numeric(Y_1))
    
    # Survey design object with weights and FPC
    combined_data <- bind_rows(ps_with_probs, overlap_with_Y1)
    # dsgn <- survey::svydesign(
    #   ids     = ~1,
    #   strata  = ~age_category,  # Stratified by age category
    #   weights = ~survey_weight,
    #   # fpc     = ~fpc_age,       # Add FPC for age-level strata
    #   # data    =  ps_with_probs
    #   data    = combined_data
    # )
    
    # Create survey design based on sampling scheme
    dsgn <- create_survey_design(combined_data, sampling_scheme)
    
    # Calculate state-level means and standard errors
    svyby(~Y_val, ~ST, dsgn, svymean, keep.var = TRUE, deff = FALSE)
  })
  
  # -----------------------------------------------------------
  # 5. Process Group-Level (State-Level) Estimates
  # -----------------------------------------------------------
  for (st in all_states) {
    # Extract true value for this state
    true_value <- true_values[[as.character(st)]]
    
    # Extract PS estimates for this state across draws
    ps_draws <- sapply(group_ps_estimates, function(res) {
      if (st %in% res$ST) res[res$ST == st, "Y_val"] else NA
    })
    ps_se <- sapply(group_ps_estimates, function(res) {
      if (st %in% res$ST) res[res$ST == st, "se"] else NA
    })
    
    # Extract MRP estimates for this state across draws
    col_idx <- which(all_states == st)
    mrp_cols <- posterior_mrp |> dplyr::select(starts_with("group_mrp_estimate"))
    nps_draws <- mrp_cols[[col_idx]]
    
    # 5a. PS Variance (Multiple Imputation Framework)
    if (sum(!is.na(ps_draws)) > 1 && sum(!is.na(ps_se)) > 1) {
      B_ps <- var(ps_draws, na.rm = TRUE)         # Between-imputation variance
      W_ps <- mean(ps_se^2, na.rm = TRUE)         # Within-imputation variance  
      T_ps <- W_ps + B_ps/num_draws               # Total variance (Reiter's adjusted formula)
      hat_ps <- mean(ps_draws, na.rm = TRUE)      # Point estimate
      
      # Calculate bias and MSE
      bias_ps <- hat_ps - true_value
      mse_ps <- T_ps + bias_ps^2
    } else {
      B_ps <- W_ps <- T_ps <- hat_ps <- bias_ps <- mse_ps <- NA
    }
    
    mse_ps <- 2
    
    # 5b. NPS/MRP Variance (Bayesian Posterior)
    hat_nps <- mean(nps_draws, na.rm = TRUE)    # Point estimate
    T_nps <- var(nps_draws, na.rm = TRUE)       # Total variance (posterior variance)
    
    # Calculate bias and MSE
    bias_nps <- hat_nps - true_value
    mse_nps <- T_nps + bias_nps^2
    
    mse_nps <- 2
    
    # 5c. Combined Estimate - DIRECT APPROACH
    if (use_mse) {
      # Use MSE-based weighting
      w_ps <- 1/mse_ps
      w_nps <- 1/mse_nps
      
      # Combined estimate based on overall means
      hat_comb <- (hat_ps * w_ps + hat_nps * w_nps) / (w_ps + w_nps)
      
      # Variance of the combined estimate - using CORRECTED MSE-based formula
      # T_comb = (T_PS * MSE_NPS^2 + T_NPS * MSE_PS^2) / (MSE_PS + MSE_NPS)^2
      T_comb <- (T_ps * mse_nps^2 + T_nps * mse_ps^2) / (mse_ps + mse_nps)^2
    } else {
      # Use precision-based weighting
      w_ps <- 1/T_ps
      w_nps <- 1/T_nps
      
      # Combined estimate based on overall means
      hat_comb <- (hat_ps * w_ps + hat_nps * w_nps) / (w_ps + w_nps)
      
      # Variance of the combined estimate - using variance-based formula
      T_comb <- (T_ps * T_nps) / (T_ps + T_nps)
    }
    
    # Calculate bias and MSE for combined estimate
    bias_comb <- hat_comb - true_value
    mse_comb <- T_comb + bias_comb^2
    
    # Store state-level results
    group_results[[as.character(st)]] <- list(
      hat_theta_ps  = hat_ps,   T_s_ps  = T_ps,   B_ps  = B_ps,   W_ps  = W_ps,
      bias_ps = bias_ps, mse_ps = mse_ps,
      
      hat_theta_nps = hat_nps,  T_s_nps = T_nps,
      bias_nps = bias_nps, mse_nps = mse_nps,
      
      hat_theta_comb = hat_comb, T_s_comb = T_comb,
      bias_comb = bias_comb, mse_comb = mse_comb
    )
  }
  
  # -----------------------------------------------------------
  # 6. Calculate Overall Estimates (Across All States)
  # -----------------------------------------------------------
  
  # Extract true value for overall
  true_value_overall <- true_values[["overall"]]
  
  # 6a. PS overall estimates across draws
  overall_ps_list <- lapply(seq_len(num_draws), function(m) {
    # Get imputed PROBABILITIES for PS-only sample
    ps_with_probs <- ps_sample |> mutate(
      Y_val = unlist(posterior_probs[m, 1:nrow(ps_sample)])
    )
    # For overlap sample, use actual observed Y1 BINARY values
    overlap_with_Y1 <- overlap_sample |> mutate(Y_val = as.numeric(Y_1))
    
    combined_data <- bind_rows(ps_with_probs, overlap_with_Y1)
    # Create survey design based on sampling scheme
    dsgn <- create_survey_design(combined_data, sampling_scheme)
    svymean(~Y_val, dsgn)
  })
  
  ps_vals <- sapply(overall_ps_list, as.numeric)
  ps_vars <- sapply(overall_ps_list, function(x) attr(x, "var"))
  
  # 6b. PS overall variance (multiple imputation framework)
  B_ps_overall <- var(ps_vals, na.rm = TRUE)             # Between-imputation variance
  W_ps_overall <- mean(ps_vars, na.rm = TRUE)            # Within-imputation variance
  T_ps_overall <- W_ps_overall + B_ps_overall/num_draws  # Total variance (Reiter's adjusted formula)
  hat_ps_overall <- mean(ps_vals, na.rm = TRUE)    # Point estimate
  
  # Calculate bias and MSE for PS overall
  bias_ps_overall <- hat_ps_overall - true_value_overall
  mse_ps_overall <- T_ps_overall + bias_ps_overall^2
  
  mse_ps_overall <- 2
  
  # 6c. NPS overall variance (Bayesian posterior)
  hat_nps_overall <- mean(overall_mrp, na.rm = TRUE)  # Point estimate
  T_nps_overall <- var(overall_mrp, na.rm = TRUE)     # Posterior variance
  
  # Calculate bias and MSE for NPS overall
  bias_nps_overall <- hat_nps_overall - true_value_overall
  mse_nps_overall <- T_nps_overall + bias_nps_overall^2
  
  mse_nps_overall <- 2
  
  # 6d. Combined overall estimate - DIRECT APPROACH
  if (use_mse) {
    # Use MSE-based weighting
    w_ps_overall <- 1/mse_ps_overall
    w_nps_overall <- 1/mse_nps_overall
    
    # Combined estimate based on overall means
    hat_comb_overall <- (hat_ps_overall * w_ps_overall + hat_nps_overall * w_nps_overall) / 
      (w_ps_overall + w_nps_overall)
    
    # Variance of the combined estimate - using CORRECTED MSE-based formula
    # T_comb = (T_PS * MSE_NPS^2 + T_NPS * MSE_PS^2) / (MSE_PS + MSE_NPS)^2
    T_comb_overall <- (T_ps_overall * mse_nps_overall^2 + T_nps_overall * mse_ps_overall^2) / 
      (mse_ps_overall + mse_nps_overall)^2
  } else {
    # Use precision-based weighting
    w_ps_overall <- 1/T_ps_overall
    w_nps_overall <- 1/T_nps_overall
    
    # Combined estimate based on overall means
    hat_comb_overall <- (hat_ps_overall * w_ps_overall + hat_nps_overall * w_nps_overall) / 
      (w_ps_overall + w_nps_overall)
    
    # Variance of the combined estimate - using variance-based formula
    T_comb_overall <- (T_ps_overall * T_nps_overall) / (T_ps_overall + T_nps_overall)
  }
  
  # Calculate bias and MSE for combined overall
  bias_comb_overall <- hat_comb_overall - true_value_overall
  mse_comb_overall <- T_comb_overall + bias_comb_overall^2
  
  # 6e. Store overall results
  group_results[["overall"]] <- list(
    hat_theta_ps   = hat_ps_overall,   T_s_ps   = T_ps_overall,   
    B_ps           = B_ps_overall,     W_ps     = W_ps_overall,
    bias_ps        = bias_ps_overall,  mse_ps    = mse_ps_overall,
    
    hat_theta_nps  = hat_nps_overall,  T_s_nps  = T_nps_overall,
    bias_nps       = bias_nps_overall, mse_nps   = mse_nps_overall,
    
    hat_theta_comb = hat_comb_overall, T_s_comb = T_comb_overall,
    bias_comb      = bias_comb_overall, mse_comb = mse_comb_overall
  )
  
  return(group_results)
}

# Function to run a full batch of simulations
# Automatically names files according to MAR or MNAR
run_OG_simulations <- function(S, pop_data, true_values, overlap_pct = 0.11, sampling_scheme,
                               use_mnar = FALSE, use_mse = TRUE) {
  # Set up parallel processing
  plan(multisession, workers = 3)  # Adjust based on available cores
  
  # Use progressr to track progress
  progressr::with_progress({
    p <- progressr::progressor(along = 1:S)
    
    # Run simulations in parallel
    simulation_results <- future_lapply(1:S, function(x) {
      p(sprintf("Simulation %d completed", x))  # Update progress bar
      run_simulation_theta1_modified(
        s = x,
        pop_data = pop_data,
        overlap_pct = overlap_pct,
        true_values = true_values,
        sampling_scheme = sampling_scheme,
        use_mnar = use_mnar,
        use_mse = use_mse
      )
    }, future.seed = TRUE)
  })
  
  # Clean up parallel workers
  plan(multisession)
  gc()
  
  # Determine MAR/MNAR suffix for file names
  mnar_suffix <- if(use_mnar) "MNAR" else "MAR"
  
  # Save results with appropriate naming
  results_filename <- paste0(sampling_scheme, "-rac3-Avg-OG-Predictive-", mnar_suffix, 
                             "-Results-", overlap_pct*100, "pct-", S, "s.rds")
  saveRDS(simulation_results, results_filename)
  cat("Results saved to:", results_filename, "\n")
  
  # Calculate metrics
  OG_metrics <- calculate_metrics_unified(
    simulation_results = simulation_results,
    true_values = true_values,
    method = c("ps")  # Make sure this matches the suffix used in group_results
  )
  
  # Save metrics with appropriate naming
  metrics_filename <- paste0(sampling_scheme, "-rac3-Avg-OG-Predictive-", mnar_suffix, 
                             "-Metrics-", overlap_pct*100, "pct-", S, "s.rds")
  saveRDS(OG_metrics, metrics_filename)
  cat("Metrics saved to:", metrics_filename, "\n")
  
  return(list(results = simulation_results, metrics = OG_metrics))
}

# Run all simulations with a loop to avoid repetition
run_all_simulations <- function(S = 100, pop_data, true_values, overlap_pct = 0.11, use_mse) {
  
  # Define simulation parameters
  sampling_schemes <- c("SRS", "Stratified")
  mnar_settings <- c(TRUE, FALSE)  # TRUE for MNAR, FALSE for MAR
  
  # Store all results
  all_results <- list()
  
  # Run all combinations
  for (scheme in sampling_schemes) {
    for (use_mnar in mnar_settings) {
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat("Running:", scheme, if(use_mnar) "MNAR" else "MAR", "\n")
      cat(paste(rep("=", 60), collapse = ""), "\n\n")
      
      result <- run_OG_simulations(
        S = S, 
        pop_data = pop_data, 
        true_values = true_values,
        use_mnar = use_mnar, 
        use_mse = use_mse,
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

# Example usage:
# Run all 6 combinations (3 schemes × 2 MAR/MNAR)
all_results <- run_all_simulations(S = 100,
                                   pop_data = pop_data_predictive,
                                   true_values = true_values_predictive,
                                   overlap_pct = 0.11,
                                   use_mse = FALSE)

# Or run individually with automatic naming:
run_OG_simulations(S = 100, pop_data = pop_data_predictive, true_values = true_values_predictive,
                   use_mnar = FALSE, use_mse = TRUE, overlap_pct = 0.11,
                   sampling_scheme = "Stratified")

calculate_metrics_unified(simulation_results = `SRS-MI-PUMS-Predictive-MAR-OG-Results-11pct-100s`,
                          true_values = true_values_predictive,
                          method = c("ps"))
calculate_metrics_unified(simulation_results = `SRS-known-N_c-MI-PUMS-Predictive-MAR-OG-Results-11pct-100s`,
                          true_values = true_values_predictive,
                          method = c("ps"))
calculate_metrics_unified(simulation_results = `SRS-MI-PUMS-Predictive-X*Y2-MAR-OG-Results-11pct-100s`,
                          true_values = true_values_predictive,
                          method = c("ps"))
calculate_metrics_unified(simulation_results = `SRS-OG-Predictive-MAR-Results-11pct-100s`,
                          true_values = true_values_predictive,
                          method = c("ps"))


calculate_metrics_unified(simulation_results = `SRS-rac3-Precision-OG-Predictive-MAR-Results-11pct-100s`,
                          true_values = true_values_predictive,
                          method = c("ps"))


calculate_metrics_unified(simulation_results = `SRS-rac3-Avg-OG-Predictive-MAR-Results-11pct-100s`,
                          true_values = true_values_predictive,
                          method = c("ps"))
