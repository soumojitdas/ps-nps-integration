library(cmdstanr)
library(rstanarm)
library(posterior)
library(tidyverse)
library(future)
library(future.apply)
library(survey)

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



samples_OG_predictive_ <- create_overlap_samples_age_only(pop_data = pop_data_predictive, 
                                                          sampling_scheme = "SRS", use_mnar = F, seed = 17,
                                                          overlap_pct = 0.11)


create_samples_MRP <- function(pop_data,
                               N_ps = 1100,
                               N_nps = 4300,
                               overlap_pct = 0.11,
                               stratify_by = "age_category",  # Default: age (what you already ran)
                               # NPS selection parameters
                               nps_betas = c(-2.703002, 0.1404789, 0.4330459, 0.4223657, 0.2616391),
                               nps_sigma = 0.1,
                               # MNAR mechanism parameters
                               use_mnar = FALSE,
                               mnar_shift = 0.5,
                               seed = 17) {
  set.seed(seed)
  
  # Validate inputs
  required_vars <- c("Y_1", "Y_2", "ST", "age_category", "race_category", "gender_binary", stratify_by)
  missing_vars <- setdiff(required_vars, names(pop_data))
  if (length(missing_vars) > 0) {
    stop("Population data missing required variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Ensure factors are properly set
  pop_data <- pop_data %>%
    mutate(
      age_category = factor(age_category, levels = c("18-29", "30-44", "45-59", "60-74", "75+")),
      race_category = factor(race_category),
      gender_binary = factor(gender_binary)
    )
  
  # Set stratification factor
  if (stratify_by %in% names(pop_data)) {
    pop_data[[stratify_by]] <- factor(pop_data[[stratify_by]])
  }
  
  cat("=== MRP SAMPLING SETUP ===\n")
  cat("Stratification variable:", stratify_by, "\n")
  cat("Model variables: age_category, race_category, gender_binary\n")
  cat("Note: MRP ignores Y2 entirely (violates DGP)\n")
  cat("MNAR mechanism:", ifelse(use_mnar, "ENABLED", "DISABLED"), "\n")
  
  # Create cells for MRP (always based on demographics, not stratification var)
  poststrat_cells <- pop_data %>%
    mutate(cell_id = interaction(ST, age_category, race_category, gender_binary, drop = TRUE)) %>%
    group_by(cell_id, ST, age_category, race_category, gender_binary) %>%
    summarise(
      n_pop = n(),
      Y_1_mean = mean(Y_1),
      Y_2_mean = mean(Y_2),
      .groups = "drop"
    ) %>%
    mutate(cell_index = as.integer(factor(cell_id)))
  
  cat("\nNumber of MRP cells (ST × age × race × gender):", nrow(poststrat_cells), "\n")
  
  # Create cell-level design matrix
  X_cell <- model.matrix(~ age_category + race_category + gender_binary, data = poststrat_cells)
  
  # Map population to cells
  pop_data <- pop_data %>%
    mutate(cell_id = interaction(ST, age_category, race_category, gender_binary, drop = TRUE)) %>%
    left_join(poststrat_cells %>% select(cell_id, cell_index, n_pop), by = "cell_id")
  
  # =====================================================
  # PROBABILITY SAMPLE - STRATIFIED SAMPLING
  # =====================================================
  
  # Get stratification variable levels
  strata <- pop_data %>%
    group_by(!!sym(stratify_by)) %>%
    summarise(N_h = n(), .groups = "drop") %>%
    mutate(n_h = round(N_ps * N_h / sum(N_h)))
  
  # Adjust to ensure exactly N_ps
  diff <- N_ps - sum(strata$n_h)
  if (diff != 0) strata$n_h[1] <- strata$n_h[1] + diff
  
  cat("\nStratified sampling allocation by", stratify_by, ":\n")
  print(strata)
  
  # Perform stratified sampling
  ps_set <- integer(0)
  for (i in 1:nrow(strata)) {
    stratum_value <- strata[[stratify_by]][i]
    target_n <- strata$n_h[i]
    stratum_indices <- which(pop_data[[stratify_by]] == stratum_value)
    
    if (target_n > length(stratum_indices)) {
      warning(sprintf("Requested %d samples from stratum %s, but only %d units available",
                      target_n, stratum_value, length(stratum_indices)))
      sampled <- stratum_indices
    } else {
      sampled <- sample(stratum_indices, target_n, replace = FALSE)
    }
    ps_set <- c(ps_set, sampled)
  }
  
  cat("PS sample size:", length(ps_set), "\n")
  
  # Weight function for stratified sampling
  assign_ps_weights <- function(data) {
    data %>%
      left_join(strata %>% select(all_of(c(stratify_by, "N_h", "n_h"))), by = stratify_by) %>%
      mutate(survey_weight = N_h / n_h) %>%
      select(-N_h, -n_h)
  }
  
  # =====================================================
  # NON-PROBABILITY SAMPLE - BIASED BY AGE
  # =====================================================
  
  cat("\n--- NPS Selection (age-biased) ---\n")
  
  # Create design matrix for NPS selection (age-based)
  X_nps_selection <- model.matrix(~ age_category, data = pop_data)
  
  # Generate errors and calculate selection probabilities
  nps_errors <- rnorm(nrow(pop_data), 0, nps_sigma)
  logit_p_nps <- X_nps_selection %*% nps_betas + nps_errors
  
  # Add MNAR effect if requested
  if (use_mnar) {
    cat("Applying MNAR mechanism based on Y_2\n")
    logit_p_nps <- logit_p_nps + mnar_shift * pop_data$Y_2
  }
  
  nps_selection_prob <- plogis(logit_p_nps)
  
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
  
  # =====================================================
  # CREATE FINAL DATASETS
  # =====================================================
  
  ps_only_set <- setdiff(ps_set, overlap)
  nps_only_set <- setdiff(nps_set, overlap)
  
  # Create datasets
  ps_data <- pop_data[ps_only_set, ] %>% mutate(Y_1 = NA)
  overlap_data <- pop_data[overlap, ]
  nps_data <- pop_data[nps_only_set, ] %>% mutate(Y_2 = NA)
  
  # Assign weights
  ps_data <- assign_ps_weights(ps_data)
  overlap_data <- assign_ps_weights(overlap_data)
  ps_with_overlap <- bind_rows(ps_data, overlap_data)
  
  # Cell statistics
  cell_summary <- ps_with_overlap %>%
    group_by(cell_index) %>%
    summarise(
      n_ps_c = n(),
      w_bar_c = mean(survey_weight),
      cell_id = first(cell_id),
      ST = first(ST),
      .groups = "drop"
    ) %>%
    left_join(poststrat_cells %>% select(cell_index, n_pop), by = "cell_index")
  
  # State mapping
  unique_states <- sort(unique(pop_data$ST))
  state_to_index <- setNames(1:length(unique_states), unique_states)
  
  # Create design matrices for MRP model
  X_nps <- model.matrix(~ age_category + race_category + gender_binary, data = nps_data)
  
  # Map cells
  nps_cell_id <- match(nps_data$cell_index, cell_summary$cell_index)
  valid_nps <- !is.na(nps_cell_id)
  
  # Create MRP-specific Stan data (simpler than OG)
  mrp_stan_data <- list(
    # NPS data for outcome model
    N_nps = sum(valid_nps),
    J = length(unique_states),
    K = ncol(X_nps),
    group_nps = state_to_index[nps_data$ST[valid_nps]],
    X_nps = X_nps[valid_nps, , drop = FALSE],
    Y1_nps = nps_data$Y_1[valid_nps],
    
    # Cell information for poststratification
    C = nrow(cell_summary),
    n_ps_c = cell_summary$n_ps_c,
    w_bar_c = cell_summary$w_bar_c,
    cell_to_group = state_to_index[cell_summary$ST],
    X_cell = X_cell[match(cell_summary$cell_index, poststrat_cells$cell_index), ],
    nps_cell_id = nps_cell_id[valid_nps]
  )
  
  # =====================================================
  # ANALYSIS OF STRATIFICATION EFFECT
  # =====================================================
  
  cat("\n=== STRATIFICATION ANALYSIS ===\n")
  
  # Show distribution of stratification variable
  strat_pop <- prop.table(table(pop_data[[stratify_by]]))
  strat_ps <- prop.table(table(ps_data[[stratify_by]]))
  strat_nps <- prop.table(table(nps_data[[stratify_by]]))
  
  cat("\nDistribution of", stratify_by, ":\n")
  cat("Population:", paste(names(strat_pop), "=", round(strat_pop, 3)), "\n")
  cat("PS sample: ", paste(names(strat_ps), "=", round(strat_ps, 3)), "\n")
  cat("NPS sample:", paste(names(strat_nps), "=", round(strat_nps, 3)), "\n")
  
  # If stratifying by income, show the effect MRP will miss
  if (stratify_by == "income_category") {
    cat("\n=== Y1 by Income (MRP cannot adjust for this!) ===\n")
    income_y1 <- pop_data %>%
      group_by(income_category) %>%
      summarise(
        Y1_mean = mean(Y_1),
        n = n(),
        prop = n() / nrow(pop_data),
        .groups = "drop"
      )
    print(income_y1)
    
    cat("\nMRP models Y1 ~ age + race + gender, completely missing income effects!\n")
  }
  
  # Return simplified structure
  return(list(
    # Data for MRP
    ps_data = ps_data,
    nps_data = nps_data[valid_nps, ],
    overlap_data = overlap_data,
    
    # Stan data
    mrp_stan_data = mrp_stan_data,
    
    # Metadata
    actual_overlap_pct = length(overlap) / length(ps_set),
    unique_states = unique_states,
    stratify_by = stratify_by,
    n_ps = length(ps_set),
    n_nps = length(nps_set),
    n_overlap = length(overlap)
  ))
}

# For comparison with your existing results (age stratification)
samples_mrp_age <- create_samples_MRP(
  pop_data = pop_data_predictive,
  stratify_by = "age_category",  # Matches what you already did
  seed = 17
)

# To show MRP failure (income stratification)
samples_mrp_income <- create_samples_MRP(
  pop_data = pop_data_predictive,  # With strong income effect
  stratify_by = "income_category",
  seed = 17
)

# MRP Predictive Model Simulation Function

# Compile the MRP predictive model once
mrp_predictive_model <- cmdstan_model("MRP_predictive.stan")

run_simulation_MRP_predictive_income <- function(s, pop_data, overlap_pct, true_values, 
                                                 stratify_by = "income_category", use_mnar = FALSE) {
  # -----------------------------------------------------------
  # 1. Create Samples and Prepare Stan Data
  # -----------------------------------------------------------
  samples <- create_samples_MRP(
    pop_data = pop_data,
    N_ps = 1100,
    N_nps = 4300,
    overlap_pct = overlap_pct,
    stratify_by = stratify_by,
    use_mnar = use_mnar,
    seed = s + 17
  )
  
  # Extract MRP-specific Stan data
  mrp_stan_data <- samples$mrp_stan_data
  
  # -----------------------------------------------------------
  # 2. Fit the MRP Stan Model
  # -----------------------------------------------------------
  fit <- mrp_predictive_model$sample(
    data = mrp_stan_data,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    seed = s + 42,
    refresh = 0
  )
  
  # -----------------------------------------------------------
  # 3. Extract Posterior Draws (Already Blended in Stan)
  # -----------------------------------------------------------
  # Extract the blended MRP estimates directly from Stan
  group_mrp_draws <- fit$draws("group_mrp_estimate", format = "df") |> as_tibble()
  overall_mrp <- fit$draws("overall_mrp_estimate", format = "draws_matrix") |> as.numeric()
  
  all_states <- sort(unique(pop_data$ST))
  group_results <- list()
  
  # -----------------------------------------------------------
  # 4. Process Group-Level (State-Level) Estimates
  # -----------------------------------------------------------
  for (st in all_states) {
    # Extract true value for this state
    true_value <- true_values[[as.character(st)]]
    
    # Get state index
    st_idx <- which(samples$unique_states == st)
    
    # Extract MRP estimates for this state across draws
    col_name <- paste0("group_mrp_estimate[", st_idx, "]")
    mrp_draws <- group_mrp_draws[[col_name]]
    
    # Filter valid draws (exclude -1 values which indicate no data)
    valid_draws <- mrp_draws[mrp_draws >= 0]
    
    if (length(valid_draws) > 0) {
      # Calculate statistics
      hat_mrp <- mean(valid_draws, na.rm = TRUE)
      T_mrp <- var(valid_draws, na.rm = TRUE)
      
      # Calculate bias and MSE
      bias_mrp <- hat_mrp - true_value
      mse_mrp <- T_mrp + bias_mrp^2
    } else {
      hat_mrp <- T_mrp <- bias_mrp <- mse_mrp <- NA
    }
    
    # Store state-level results
    group_results[[as.character(st)]] <- list(
      hat_theta_mrp = hat_mrp,
      T_s_mrp = T_mrp,
      bias_mrp = bias_mrp,
      mse_mrp = mse_mrp
    )
  }
  
  # -----------------------------------------------------------
  # 5. Process Overall Estimate
  # -----------------------------------------------------------
  true_value_overall <- true_values[["overall"]]
  
  # Filter valid draws
  valid_overall <- overall_mrp[overall_mrp >= 0]
  
  # Calculate statistics
  hat_mrp_overall <- mean(valid_overall, na.rm = TRUE)
  T_mrp_overall <- var(valid_overall, na.rm = TRUE)
  
  # Calculate bias and MSE
  bias_mrp_overall <- hat_mrp_overall - true_value_overall
  mse_mrp_overall <- T_mrp_overall + bias_mrp_overall^2
  
  # Store overall results
  group_results[["overall"]] <- list(
    hat_theta_mrp = hat_mrp_overall,
    T_s_mrp = T_mrp_overall,
    bias_mrp = bias_mrp_overall,
    mse_mrp = mse_mrp_overall
  )
  
  return(group_results)
}

# Function to run all MRP predictive simulations with income stratification
run_MRP_predictive_income_simulations <- function(S, pop_data, true_values, overlap_pct = 0.11, 
                                                  stratify_by = "income_category", use_mnar = FALSE) {
  # Compile the Stan model once
  cat("Compiling MRP predictive Stan model...\n")
  mrp_predictive_model <<- cmdstan_model("MRP_predictive.stan")
  
  # Set up parallel processing
  plan(multisession, workers = 3)  # Adjust based on available cores
  
  # Use progressr to track progress
  progressr::with_progress({
    p <- progressr::progressor(along = 1:S)
    
    # Run simulations in parallel
    simulation_results <- future_lapply(1:S, function(x) {
      p(sprintf("Simulation %d completed", x))  # Update progress bar
      run_simulation_MRP_predictive_income(
        s = x,
        pop_data = pop_data,
        overlap_pct = overlap_pct,
        true_values = true_values,
        stratify_by = stratify_by,
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
  results_filename <- paste0(stratify_by, "-Stratified-MRP-Predictive-", mnar_suffix, 
                             "-Results-", overlap_pct*100, "pct-", S, "s.rds")
  saveRDS(simulation_results, results_filename)
  cat("Results saved to:", results_filename, "\n")
  
  # Calculate metrics using existing function
  mrp_metrics <- calculate_metrics_unified(
    simulation_results = simulation_results,
    true_values = true_values,
    method = "mrp"
  )
  
  # Save metrics with appropriate naming
  metrics_filename <- paste0(stratify_by, "-Stratified-MRP-Predictive-", mnar_suffix, 
                             "-Metrics-", overlap_pct*100, "pct-", S, "s.rds")
  saveRDS(mrp_metrics, metrics_filename)
  cat("Metrics saved to:", metrics_filename, "\n")
  
  return(list(results = simulation_results, metrics = mrp_metrics))
}

run_MRP_predictive_income_simulations(S = 100, 
                                      pop_data = pop_data_predictive_mrp,
                                      true_values = true_values_predictive_mrp,
                                      use_mnar = T)
