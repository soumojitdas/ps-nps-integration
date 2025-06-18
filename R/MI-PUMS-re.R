# Download and Process Michigan PUMS Data for Small Area Estimation ----
# This script downloads PUMS data and creates categorical variables matching your setup

library(tidycensus)
library(tidyverse)

# Set your Census API key (replace with your actual key)
# census_api_key("YOUR_API_KEY_HERE", install = TRUE)

# ========================================
# STEP 1: Download Michigan PUMS Data
# ========================================

message("Downloading Michigan PUMS data from 2021 5-year ACS...")

mi_pums_raw <- get_pums(
  variables = c(
    "PUMA",      # Geographic identifier
    "AGEP",      # Age
    "SEX",       # Gender  
    "RAC1P",     # Race (detailed)
    "HISP",      # Hispanic origin
    "SCHL",      # Educational attainment
    "PINCP",     # Personal income
    "PWGTP"      # Person weight
  ),
  state = "MI",
  year = 2021,    # Using 2021 to avoid PUMA boundary issues
  survey = "acs5",
  rep_weights = "person",  # Include replicate weights
  recode = TRUE   # Get human-readable labels
)

message(paste("Downloaded", nrow(mi_pums_raw), "person records"))

# ========================================
# STEP 2: Create Categorical Variables
# ========================================

message("Creating categorical variables...")

mi_pums_clean <- mi_pums_raw %>%
  # Filter to adults only (18+)
  filter(AGEP >= 18) %>%
  
  # Create age categories matching your setup
  mutate(
    age_category = case_when(
      AGEP >= 18 & AGEP <= 29 ~ "18-29",
      AGEP >= 30 & AGEP <= 44 ~ "30-44",
      AGEP >= 45 & AGEP <= 59 ~ "45-59",
      AGEP >= 60 & AGEP <= 74 ~ "60-74",
      AGEP >= 75 ~ "75+",
      TRUE ~ NA_character_
    ),
    age_category = factor(age_category, 
                          levels = c("18-29", "30-44", "45-59", "60-74", "75+"))
  ) %>%
  
  # Create binary gender variable
  mutate(
    gender_binary = case_when(
      SEX == 1 ~ "Male",
      SEX == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    gender_binary = factor(gender_binary, levels = c("Male", "Female"))
  ) %>%
  
  # Create race/ethnicity categories
  # Note: In PUMS, HISP == 1 means "Not Spanish/Hispanic/Latino"
  mutate(
    race_category = case_when(
      HISP > 1 ~ "Hispanic",   # HISP values 2-24 are Hispanic origins
      RAC1P == 1 ~ "White",    # White alone, non-Hispanic
      RAC1P == 2 ~ "Black",    # Black alone, non-Hispanic
      RAC1P == 6 ~ "Asian",    # Asian alone, non-Hispanic
      TRUE ~ "Other"           # All other races, non-Hispanic
    ),
    race_category = factor(race_category, 
                           levels = c("White", "Black", "Hispanic", "Asian", "Other"))
  ) %>%
  
  # Create education categories
  mutate(
    education = case_when(
      SCHL <= 15 ~ "Less than HS",           # No schooling to Grade 11
      SCHL %in% 16:17 ~ "High school",       # HS diploma or GED
      SCHL %in% 18:20 ~ "Some college",      # Some college or Associate's
      SCHL >= 21 ~ "Bachelor's or higher",   # Bachelor's or higher
      TRUE ~ NA_character_
    ),
    education = factor(education,
                       levels = c("Less than HS", "High school", 
                                  "Some college", "Bachelor's or higher"))
  ) %>%
  
  # Create income categories
  mutate(
    # Handle missing income (negative values in PUMS)
    income_clean = ifelse(PINCP < 0, NA, PINCP),
    income_category = case_when(
      is.na(income_clean) ~ "Missing",
      income_clean < 35000 ~ "<$35k",
      income_clean < 55000 ~ "$35-55k",
      income_clean < 100000 ~ "$55-100k",
      income_clean >= 100000 ~ ">$100k",
      TRUE ~ NA_character_
    ),
    income_category = factor(income_category,
                             levels = c("<$35k", "$35-55k", "$55-100k", ">$100k", "Missing"))
  )

# ========================================
# STEP 3: Select 10 PUMAs for Analysis
# ========================================

# Get PUMA population sizes
puma_sizes <- mi_pums_clean %>%
  group_by(PUMA) %>%
  summarise(
    n_sample = n(),
    pop_estimate = sum(PWGTP),
    .groups = "drop"
  ) %>%
  arrange(desc(pop_estimate))

message("\nTop 20 PUMAs by population:")
print(puma_sizes %>% head(20))

# Select 10 diverse PUMAs
# Strategy: Mix of large (urban) and medium-sized PUMAs
selected_pumas <- c(
  "03211",  # Detroit City (South Central & Southeast)
  "03209",  # Detroit City (North Central)
  "01703",  # Flint City Area
  "02702",  # Ann Arbor City Area
  "01802",  # Ingham County (Northwest) - Lansing
  "01002",  # Grand Rapids City Area
  "02102",  # Kalamazoo City Area
  "02906",  # Birmingham & Bloomfield Area
  "03003",  # Sterling Heights
  "00100"   # Western Upper Peninsula
)

# Filter to selected PUMAs
mi_pums_selected <- mi_pums_clean %>%
  filter(PUMA %in% selected_pumas)

message(paste("\nFiltered to", n_distinct(mi_pums_selected$PUMA), 
              "PUMAs with", nrow(mi_pums_selected), "person records"))

# ========================================
# STEP 4: Create Analysis Dataset
# ========================================

# Create final dataset with renamed variables matching your setup
mi_analysis_data <- mi_pums_selected %>%
  transmute(
    # Individual ID
    person_id = row_number(),
    
    # Geographic identifier (PUMA as "state" in your notation)
    ST = paste0("PUMA_", PUMA),
    
    # Demographics (matching your X matrix)
    age_category = age_category,
    race_category = race_category,
    gender_binary = gender_binary,
    
    # Additional variables
    education = education,
    income_category = income_category,
    
    # Continuous variables (kept for reference)
    age_continuous = AGEP,
    income_continuous = income_clean,
    
    # Survey weight
    survey_weight = PWGTP,
    
    # Keep PUMA for reference
    PUMA_code = PUMA
  ) %>%
  # Remove any records with missing key demographics
  filter(
    !is.na(age_category),
    !is.na(race_category),
    !is.na(gender_binary)
  )

# ========================================
# STEP 5: Summary Statistics
# ========================================

message("\n=== SUMMARY STATISTICS ===")

# Overall sample size
message(paste("Total sample size:", nrow(mi_analysis_data)))
message(paste("Population represented:", 
              format(sum(mi_analysis_data$survey_weight), big.mark = ",")))

# By PUMA
message("\nSample sizes by PUMA:")
mi_analysis_data %>%
  group_by(ST, PUMA_code) %>%
  summarise(
    n_sample = n(),
    pop_size = sum(survey_weight),
    .groups = "drop"
  ) %>%
  print()

# Demographics distribution
message("\nAge distribution:")
print(table(mi_analysis_data$age_category))

message("\nRace distribution:")
print(table(mi_analysis_data$race_category))

message("\nGender distribution:")
print(table(mi_analysis_data$gender_binary))

# ========================================
# STEP 6: Save Processed Data
# ========================================

# Save the analysis dataset
saveRDS(mi_analysis_data, "mi_pums_analysis_data.rds")
message("\nSaved processed data to 'mi_pums_analysis_data.rds'")

# Create a population summary for reference
pop_summary <- mi_analysis_data %>%
  group_by(ST) %>%
  summarise(
    n = n(),
    pop_size = sum(survey_weight),
    mean_age = weighted.mean(age_continuous, survey_weight),
    pct_female = weighted.mean(gender_binary == "Female", survey_weight),
    pct_white = weighted.mean(race_category == "White", survey_weight),
    .groups = "drop"
  )

print(pop_summary)

# Generate Y_1 and Y_2 outcomes for Michigan PUMS data ----
# Following the current DGP structure from gen-pop-Li.R

mi_analysis_data <- readRDS("mi_pums_analysis_data.rds")
glimpse(mi_analysis_data)

table(mi_analysis_data$race_category)
table(mi_analysis_data$gender_binary)
table(mi_analysis_data$age_category)

# Create a new collapsed race category
mi_analysis_data <- mi_analysis_data %>%
  mutate(race_category = case_when(
    race_category == "White" ~ "White",
    race_category == "Black" ~ "Black", 
    race_category %in% c("Hispanic", "Asian", "Other") ~ "Others",
    TRUE ~ NA_character_  # Handle any unexpected values
  ))

# Convert to factor if desired
mi_analysis_data$race_category <- factor(
  mi_analysis_data$race_category,
  levels = c("White", "Black", "Others")
)

generate_outcomes_michigan <- function(mi_data, seed = 123) {
  set.seed(seed)
  
  # Helper functions
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  # Get number of observations and states (PUMAs)
  n_population <- nrow(mi_data)
  unique_states <- unique(mi_data$ST)
  n_states <- length(unique_states)
  
  cat("Generating outcomes for", n_population, "individuals in", n_states, "PUMAs\n")
  
  # ========================================
  # PARAMETERS (matching gen-pop-Li.R structure)
  # ========================================
  
  # Small area (PUMA) random effects
  sigma_u <- 0.15  # SD of random effects
  u_j <- rnorm(n_states, 0, sigma_u)
  names(u_j) <- unique_states
  
  # Prevalence model parameters (for Y_1)
  beta_intercept <- -1.5
  
  # Using the same parameter values from gen-pop-Li.R
  # but adjusting for actual Michigan demographics
  beta_race <- c(0.00, 1.9071626, 1.1448769, -0.7645307, -1.4574325)  # White (ref), Black, Hispanic, Asian, Other
  beta_gender <- c(0.00, 0.24)  # Male (ref), Female
  beta_age <- c(0.00, -1.093468881, 0.295241218, 0.006885942, 1.157410886)  # 18-29 (ref), 30-44, 45-59, 60-74, 75+
  
  # Measurement error parameters
  alpha_sens <- 1.5  # Sensitivity on logit scale
  alpha_spec <- 2.0  # Specificity on logit scale
  sens_prob <- inv_logit(alpha_sens)  # ~0.818
  spec_prob <- inv_logit(alpha_spec)  # ~0.881
  
  cat("\nModel parameters:\n")
  cat("Sensitivity:", round(sens_prob, 3), "\n")
  cat("Specificity:", round(spec_prob, 3), "\n")
  cat("Random effects SD:", sigma_u, "\n")
  
  # ========================================
  # CREATE DESIGN MATRIX
  # ========================================
  
  # Ensure factors are properly set
  mi_data <- mi_data %>%
    mutate(
      race_category = factor(race_category, 
                             levels = c("White", "Black", "Hispanic", "Asian", "Other")),
      gender_binary = factor(gender_binary, 
                             levels = c("Male", "Female")),
      age_category = factor(age_category,
                            levels = c("18-29", "30-44", "45-59", "60-74", "75+"))
    )
  
  # Create design matrix (same order as Stan model expects)
  X <- model.matrix(~ race_category + gender_binary + age_category, data = mi_data)
  
  # Combine all beta coefficients in the correct order
  beta_all <- c(beta_intercept, beta_race[-1], beta_gender[-1], beta_age[-1])
  
  # ========================================
  # GENERATE Y_1 (True Outcome)
  # ========================================
  
  cat("\nGenerating Y_1 (true outcome)...\n")
  
  # Calculate linear predictor
  eta_prev <- as.vector(X %*% beta_all)
  
  
  # Add PUMA random effects
  for (i in 1:n_population) {
    puma <- mi_data$ST[i]
    eta_prev[i] <- eta_prev[i] + u_j[puma]
  }
  
  # Generate Y_1
  prob_Y1 <- inv_logit(eta_prev)
  Y_1 <- rbinom(n_population, 1, prob_Y1)
  
  # ========================================
  # GENERATE Y_2 (Proxy Measure)
  # ========================================
  
  cat("Generating Y_2 (proxy measure)...\n")
  
  # Y_2 depends on Y_1 through measurement error model
  Y_2 <- integer(n_population)
  
  for (i in 1:n_population) {
    if (Y_1[i] == 1) {
      # True positive rate = sensitivity
      Y_2[i] <- rbinom(1, 1, sens_prob)
    } else {
      # False positive rate = 1 - specificity
      Y_2[i] <- rbinom(1, 1, 1 - spec_prob)
    }
  }
  
  # ========================================
  # ADD TO DATASET
  # ========================================
  
  mi_data_with_outcomes <- mi_data %>%
    mutate(
      Y_1 = Y_1,
      Y_2 = Y_2
    )
  
  # ========================================
  # SUMMARY STATISTICS
  # ========================================
  
  cat("\n=== OUTCOME GENERATION SUMMARY ===\n")
  
  # Overall prevalences
  cat("\nOverall prevalences:\n")
  cat("Y_1 (true outcome):", round(mean(Y_1), 3), "\n")
  cat("Y_2 (proxy measure):", round(mean(Y_2), 3), "\n")
  cat("Correlation(Y_1, Y_2):", round(cor(Y_1, Y_2), 3), "\n")
  
  # By PUMA
  cat("\nPrevalence by PUMA:\n")
  puma_summary <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(
      n = n(),
      Y_1_prev = mean(Y_1),
      Y_2_prev = mean(Y_2),
      .groups = "drop"
    ) %>%
    arrange(ST)
  print(puma_summary, n = 10)
  
  # By demographics
  cat("\nY_1 prevalence by age:\n")
  age_summary <- mi_data_with_outcomes %>%
    group_by(age_category) %>%
    summarise(
      n = n(),
      Y_1_prev = mean(Y_1),
      Y_2_prev = mean(Y_2),
      .groups = "drop"
    )
  print(age_summary)
  
  cat("\nY_1 prevalence by race:\n")
  race_summary <- mi_data_with_outcomes %>%
    group_by(race_category) %>%
    summarise(
      n = n(),
      Y_1_prev = mean(Y_1),
      Y_2_prev = mean(Y_2),
      .groups = "drop"
    )
  print(race_summary)
  
  # Measurement error validation
  cat("\n=== MEASUREMENT ERROR VALIDATION ===\n")
  
  # Calculate empirical sensitivity and specificity
  emp_sens <- sum(Y_1 == 1 & Y_2 == 1) / sum(Y_1 == 1)
  emp_spec <- sum(Y_1 == 0 & Y_2 == 0) / sum(Y_1 == 0)
  
  cat("Empirical sensitivity:", round(emp_sens, 3), 
      "(target:", round(sens_prob, 3), ")\n")
  cat("Empirical specificity:", round(emp_spec, 3), 
      "(target:", round(spec_prob, 3), ")\n")
  
  # Calculate true values for evaluation
  true_values <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(mean_Y1 = mean(Y_1), .groups = "drop") %>%
    deframe()
  
  # Add overall
  true_values <- c(true_values, overall = mean(Y_1))
  
  # Return results
  return(list(
    data = mi_data_with_outcomes,
    true_values = true_values,
    parameters = list(
      # beta = beta_all,
      beta = beta_age,
      u_j = u_j,
      sigma_u = sigma_u,
      alpha_sens = alpha_sens,
      alpha_spec = alpha_spec,
      sens_prob = sens_prob,
      spec_prob = spec_prob
    )
  ))
}

mi_data_complete <- generate_outcomes_michigan(mi_analysis_data, seed = 2025)

# Extract the data with outcomes
pop_data <- mi_data_complete$data
true_values <- mi_data_complete$true_values
true_params <- mi_data_complete$parameters

# AGE ONLY ----

generate_outcomes_michigan_age_only <- function(mi_data, seed = 123) {
  set.seed(seed)
  
  # Helper functions
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  # Get number of observations and states (PUMAs)
  n_population <- nrow(mi_data)
  unique_states <- unique(mi_data$ST)
  n_states <- length(unique_states)
  
  cat("Generating outcomes for", n_population, "individuals in", n_states, "PUMAs\n")
  
  # ========================================
  # PARAMETERS FOR AGE-ONLY MODEL
  # ========================================
  
  # Small area (PUMA) random effects
  sigma_u <- 0.15  # SD of random effects
  u_j <- rnorm(n_states, 0, sigma_u)
  names(u_j) <- unique_states
  
  # Prevalence model parameters (for Y_1) - AGE ONLY
  beta_intercept <- -1.5
  
  # Age coefficients (matching Stan expectation)
  # 18-29 is reference (implicit 0)
  beta_age <- c(0.00, -1.093, 0.295, 0.007, 1.157)  # 18-29 (ref), 30-44, 45-59, 60-74, 75+
  
  # Measurement error parameters (same as before)
  alpha_sens <- 1.5  # Sensitivity on logit scale
  alpha_spec <- 2.0  # Specificity on logit scale
  sens_prob <- inv_logit(alpha_sens)  # ~0.818
  spec_prob <- inv_logit(alpha_spec)  # ~0.881
  
  cat("\nModel parameters:\n")
  cat("Sensitivity:", round(sens_prob, 3), "\n")
  cat("Specificity:", round(spec_prob, 3), "\n")
  cat("Random effects SD:", sigma_u, "\n")
  cat("Beta (intercept):", beta_intercept, "\n")
  cat("Beta (age effects):", beta_age[-1], "\n")  # Excluding reference
  
  # ========================================
  # CREATE DESIGN MATRIX - AGE ONLY
  # ========================================
  
  # Ensure age factor is properly set
  mi_data <- mi_data %>%
    mutate(
      age_category = factor(age_category,
                            levels = c("18-29", "30-44", "45-59", "60-74", "75+"))
    )
  
  # Create design matrix - ONLY AGE (matching Stan expectation)
  X <- model.matrix(~ age_category, data = mi_data)
  
  # Check dimensions
  cat("\nDesign matrix dimensions:", dim(X), "\n")
  cat("Column names:", colnames(X), "\n")
  
  # Combine coefficients: intercept + age effects (excluding reference)
  beta_full <- c(beta_intercept, beta_age[-1])
  
  # Verify dimensions match
  if (ncol(X) != length(beta_full)) {
    stop("Design matrix columns (", ncol(X), ") don't match beta length (", 
         length(beta_full), ")")
  }
  
  # ========================================
  # GENERATE Y_1 (True Outcome)
  # ========================================
  
  cat("\nGenerating Y_1 (true outcome)...\n")
  
  # Calculate linear predictor
  eta_prev <- as.vector(X %*% beta_full)
  
  # Add PUMA random effects
  for (i in 1:n_population) {
    puma <- mi_data$ST[i]
    eta_prev[i] <- eta_prev[i] + u_j[puma]
  }
  
  # Generate Y_1
  prob_Y1 <- inv_logit(eta_prev)
  Y_1 <- rbinom(n_population, 1, prob_Y1)
  
  # ========================================
  # GENERATE Y_2 (Proxy Measure)
  # ========================================
  
  cat("Generating Y_2 (proxy measure)...\n")
  
  # Y_2 depends on Y_1 through measurement error model
  Y_2 <- integer(n_population)
  
  for (i in 1:n_population) {
    if (Y_1[i] == 1) {
      # True positive rate = sensitivity
      Y_2[i] <- rbinom(1, 1, sens_prob)
    } else {
      # False positive rate = 1 - specificity
      Y_2[i] <- rbinom(1, 1, 1 - spec_prob)
    }
  }
  
  # ========================================
  # ADD TO DATASET
  # ========================================
  
  mi_data_with_outcomes <- mi_data %>%
    mutate(
      Y_1 = Y_1,
      Y_2 = Y_2
    )
  
  # ========================================
  # SUMMARY STATISTICS
  # ========================================
  
  cat("\n=== OUTCOME GENERATION SUMMARY ===\n")
  
  # Overall prevalences
  cat("\nOverall prevalences:\n")
  cat("Y_1 (true outcome):", round(mean(Y_1), 3), "\n")
  cat("Y_2 (proxy measure):", round(mean(Y_2), 3), "\n")
  cat("Correlation(Y_1, Y_2):", round(cor(Y_1, Y_2), 3), "\n")
  
  # By PUMA (showing top 10)
  cat("\nPrevalence by PUMA (top 10):\n")
  puma_summary <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(
      n = n(),
      Y_1_prev = mean(Y_1),
      Y_2_prev = mean(Y_2),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
  print(puma_summary, n = 10)
  
  # By age (the only covariate in the model)
  cat("\nY_1 prevalence by age:\n")
  age_summary <- mi_data_with_outcomes %>%
    group_by(age_category) %>%
    summarise(
      n = n(),
      Y_1_prev = mean(Y_1),
      Y_2_prev = mean(Y_2),
      .groups = "drop"
    )
  print(age_summary)
  
  # Measurement error validation
  cat("\n=== MEASUREMENT ERROR VALIDATION ===\n")
  
  # Calculate empirical sensitivity and specificity
  emp_sens <- sum(Y_1 == 1 & Y_2 == 1) / sum(Y_1 == 1)
  emp_spec <- sum(Y_1 == 0 & Y_2 == 0) / sum(Y_1 == 0)
  
  cat("Empirical sensitivity:", round(emp_sens, 3), 
      "(target:", round(sens_prob, 3), ")\n")
  cat("Empirical specificity:", round(emp_spec, 3), 
      "(target:", round(spec_prob, 3), ")\n")
  
  # ========================================
  # CALCULATE TRUE VALUES FOR EVALUATION
  # ========================================
  
  # True values by PUMA
  true_values <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(mean_Y1 = mean(Y_1), .groups = "drop") %>%
    deframe()
  
  # Add overall
  true_values <- c(true_values, overall = mean(Y_1))
  
  # ========================================
  # RETURN RESULTS
  # ========================================
  
  return(list(
    data = mi_data_with_outcomes,
    true_values = true_values,
    parameters = list(
      beta_intercept = beta_intercept,
      beta_age = beta_age,
      beta_full = beta_full,  # What actually gets used: intercept + age[-1]
      u_j = u_j,
      sigma_u = sigma_u,
      alpha_sens = alpha_sens,
      alpha_spec = alpha_spec,
      sens_prob = sens_prob,
      spec_prob = spec_prob
    ),
    model_info = list(
      design_matrix_cols = ncol(X),
      design_matrix_names = colnames(X),
      n_population = n_population,
      n_pumas = n_states
    )
  ))
}

# ========================================
# USAGE EXAMPLE --- AGE ONLY
# ========================================

# Load your Michigan data
mi_analysis_data <- readRDS("mi_pums_analysis_data.rds")

# Generate outcomes with age-only model
mi_data_age <- generate_outcomes_michigan_age_only(mi_analysis_data, seed = 2025)

# Extract the data with outcomes
pop_data_age <- mi_data_age$data
true_values_age <- mi_data_age$true_values
true_params_age <- mi_data_age$parameters


# Now ready to use with create_overlap_samples_simplified
saveRDS(pop_data, "mi_population_with_outcomes.rds")
saveRDS(true_values, "mi_true_values.rds")
saveRDS(true_params, "mi_true_parameters.rds")


# New DGP: Y2 as Predictor Model with Random Parameter Generation ----
# Y2 is generated first as an auxiliary variable
# Y1 depends on both demographics and Y2

library(tidyverse)

generate_predictive_outcomes <- function(mi_data, seed = 123) {
  set.seed(seed)
  
  # Helper functions
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  # Get basic info
  n_population <- nrow(mi_data)
  unique_states <- sort(unique(mi_data$ST))
  n_states <- length(unique_states)
  
  cat("=== GENERATING Y2 AS PREDICTOR MODEL ===\n")
  cat("Population size:", n_population, "\n")
  cat("Number of PUMAs:", n_states, "\n")
  
  # # Ensure factors are properly set
  # mi_data <- mi_data %>%
  #   mutate(
  #     race_category = factor(race_category, 
  #                            levels = c("White", "Black", "Hispanic", "Asian", "Other")),
  #     gender_binary = factor(gender_binary, 
  #                            levels = c("Male", "Female")),
  #     age_category = factor(age_category,
  #                           levels = c("18-29", "30-44", "45-59", "60-74", "75+"))
  #   )
  
  # Get number of categories dynamically
  race_categories <- levels(factor(mi_data$race_category))
  gender_categories <- levels(factor(mi_data$gender_binary))
  age_categories <- levels(mi_data$age_category)
  
  n_race_categories <- length(race_categories)
  n_gender_categories <- length(gender_categories)
  n_age_categories <- length(age_categories)
  
  cat("\nCategory counts:\n")
  cat("- Race categories:", n_race_categories, "\n")
  cat("- Gender categories:", n_gender_categories, "\n") 
  cat("- Age categories:", n_age_categories, "\n")
  
  # Design matrix - using full demographics
  X <- model.matrix(~ race_category + gender_binary + age_category, data = mi_data)
  cat("\nDesign matrix dimensions:", dim(X), "\n")
  cat("Column names:", colnames(X), "\n")
  
  # ========================================
  # STEP 1: GENERATE Y2 (Auxiliary Variable)
  # ========================================
  
  cat("\n--- Generating Y2 (auxiliary variable) ---\n")
  
  # Parameters for Y2 model - RANDOM DRAWS
  gamma_intercept <- -1.0
  gamma_race <- c(0, rnorm(n_race_categories - 1, mean = 0, sd = 0.5))
  gamma_gender <- c(0, rnorm(n_gender_categories - 1, mean = 0, sd = 0.3))
  gamma_age <- c(0, rnorm(n_age_categories - 1, mean = 0, sd = 0.4))
  
  # Print drawn parameters
  cat("\nDrawn Y2 parameters:\n")
  cat("gamma_race:", round(gamma_race, 3), "\n")
  cat("gamma_gender:", round(gamma_gender, 3), "\n")
  cat("gamma_age:", round(gamma_age, 3), "\n")
  
  # Combine gamma coefficients in correct order matching design matrix
  gamma_all <- c(gamma_intercept, gamma_race[-1], gamma_gender[-1], gamma_age[-1])
  names(gamma_all) <- colnames(X)
  
  # Random effects for Y2
  sigma_v <- 0.2
  v_j <- rnorm(n_states, 0, sigma_v)
  names(v_j) <- unique_states
  
  # Linear predictor for Y2
  eta_Y2 <- as.vector(X %*% gamma_all)
  
  # Add state random effects
  for (i in 1:n_population) {
    state <- mi_data$ST[i]
    eta_Y2[i] <- eta_Y2[i] + v_j[state]
  }
  
  # Generate Y2
  prob_Y2 <- inv_logit(eta_Y2)
  Y_2 <- rbinom(n_population, 1, prob_Y2)
  
  cat("\nY2 prevalence:", round(mean(Y_2), 3), "\n")
  
  # ========================================
  # STEP 2: GENERATE Y1 (Outcome depending on Y2)
  # ========================================
  
  cat("\n--- Generating Y1 (outcome) ---\n")
  
  # Parameters for Y1 model - RANDOM DRAWS (independent from Y2 parameters)
  beta_intercept <- -1.5
  beta_race <- c(0, rnorm(n_race_categories - 1, mean = 0, sd = 0.6))
  beta_gender <- c(0, rnorm(n_gender_categories - 1, mean = 0, sd = 0.3))
  beta_age <- c(0, rnorm(n_age_categories - 1, mean = 0, sd = 0.5))
  
  # KEY PARAMETER: Effect of Y2 on Y1 (ensure positive)
  alpha <- abs(rnorm(1, mean = 1.5, sd = 0.3))
  
  # Print drawn parameters
  cat("\nDrawn Y1 parameters:\n")
  cat("beta_race:", round(beta_race, 3), "\n")
  cat("beta_gender:", round(beta_gender, 3), "\n")
  cat("beta_age:", round(beta_age, 3), "\n")
  cat("alpha (Y2->Y1 effect):", round(alpha, 3), "\n")
  
  # Combine beta coefficients
  beta_all <- c(beta_intercept, beta_race[-1], beta_gender[-1], beta_age[-1])
  names(beta_all) <- colnames(X)
  
  # Random effects for Y1
  sigma_u <- 0.15
  u_j <- rnorm(n_states, 0, sigma_u)
  names(u_j) <- unique_states
  
  # Create extended design matrix including Y2
  X_extended <- cbind(X, Y2 = Y_2)
  beta_extended <- c(beta_all, alpha)
  
  # Linear predictor for Y1
  eta_Y1 <- as.vector(X_extended %*% beta_extended)
  
  # Add state random effects
  for (i in 1:n_population) {
    state <- mi_data$ST[i]
    eta_Y1[i] <- eta_Y1[i] + u_j[state]
  }
  
  # Generate Y1
  prob_Y1 <- inv_logit(eta_Y1)
  Y_1 <- rbinom(n_population, 1, prob_Y1)
  
  cat("Y1 prevalence:", round(mean(Y_1), 3), "\n")
  
  # ========================================
  # STEP 3: ANALYZE RELATIONSHIP
  # ========================================
  
  cat("\n--- Relationship Analysis ---\n")
  
  # Cross-tabulation
  cross_tab <- table(Y_1, Y_2)
  cat("\nCross-tabulation of Y1 and Y2:\n")
  print(cross_tab)
  
  # Calculate implied "sensitivity" and "specificity"
  # These are NOT measurement error parameters, but describe the association
  implied_sens <- sum(Y_1 == 1 & Y_2 == 1) / sum(Y_1 == 1)
  implied_spec <- sum(Y_1 == 0 & Y_2 == 0) / sum(Y_1 == 0)
  
  cat("\nImplied association measures:\n")
  cat("P(Y2=1|Y1=1):", round(implied_sens, 3), "\n")
  cat("P(Y2=0|Y1=0):", round(implied_spec, 3), "\n")
  cat("Correlation(Y1,Y2):", round(cor(Y_1, Y_2), 3), "\n")
  
  # By demographics
  cat("\nY1 and Y2 prevalence by demographics:\n")
  
  # By age
  age_summary <- mi_data %>%
    mutate(Y_1 = Y_1, Y_2 = Y_2) %>%
    group_by(age_category) %>%
    summarise(
      n = n(),
      Y1_prev = mean(Y_1),
      Y2_prev = mean(Y_2),
      .groups = "drop"
    )
  cat("\nBy age:\n")
  print(age_summary)
  
  # By race
  race_summary <- mi_data %>%
    mutate(Y_1 = Y_1, Y_2 = Y_2) %>%
    group_by(race_category) %>%
    summarise(
      n = n(),
      Y1_prev = mean(Y_1),
      Y2_prev = mean(Y_2),
      .groups = "drop"
    )
  cat("\nBy race:\n")
  print(race_summary)
  
  # ========================================
  # STEP 4: PREPARE OUTPUT
  # ========================================
  
  # Add outcomes to data
  mi_data_with_outcomes <- mi_data %>%
    mutate(
      Y_1 = Y_1,
      Y_2 = Y_2
    )
  
  # Calculate true values by state
  true_values <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(mean_Y1 = mean(Y_1), .groups = "drop") %>%
    deframe()
  
  # Add overall
  true_values <- c(true_values, overall = mean(Y_1))
  
  # Store all parameters
  true_params <- list(
    # Y2 model parameters
    gamma = gamma_all,
    gamma_intercept = gamma_intercept,
    gamma_race = gamma_race,
    gamma_gender = gamma_gender,
    gamma_age = gamma_age,
    v_j = v_j,
    sigma_v = sigma_v,
    
    # Y1 model parameters
    beta = beta_all,
    beta_intercept = beta_intercept,
    beta_race = beta_race,
    beta_gender = beta_gender,
    beta_age = beta_age,
    alpha = alpha,  # KEY: Y2 -> Y1 effect
    u_j = u_j,
    sigma_u = sigma_u,
    
    # For comparison with measurement error interpretation
    implied_sens = implied_sens,
    implied_spec = implied_spec,
    
    # Store category information
    n_race_categories = n_race_categories,
    n_gender_categories = n_gender_categories,
    n_age_categories = n_age_categories,
    race_categories = race_categories,
    gender_categories = gender_categories,
    age_categories = age_categories
  )
  
  cat("\n=== GENERATION COMPLETE ===\n")
  cat("Key parameter α (Y2->Y1 effect):", round(alpha, 3), "\n")
  cat("This is fundamentally different from measurement error!\n")
  
  return(list(
    data = mi_data_with_outcomes,
    true_values = true_values,
    true_params = true_params,
    X = X
  ))
}

# Example usage
# Assuming mi_analysis_data is loaded
predictive_result <- generate_predictive_outcomes(mi_analysis_data, seed = 2024)

pop_data_predictive <- predictive_result$data
true_values_predictive <- predictive_result$true_values
true_params_predictive <- predictive_result$true_params

# Save for use with new Stan model
saveRDS(pop_data_predictive, "mi_population_predictive_model.rds")
saveRDS(true_values_predictive, "mi_true_values_predictive.rds")
saveRDS(true_params_predictive, "mi_true_params_predictive.rds")


# Simplified DGP: Y2 as Predictor Model with Only Age Variable
# Y2 is generated first as an auxiliary variable
# Y1 depends on age and Y2

generate_predictive_outcomes_age_only <- function(mi_data, seed = 123) {
  set.seed(seed)
  
  # Helper functions
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  # Get basic info
  n_population <- nrow(mi_data)
  unique_states <- sort(unique(mi_data$ST))
  n_states <- length(unique_states)
  
  cat("=== GENERATING Y2 AS PREDICTOR MODEL (AGE ONLY) ===\n")
  cat("Population size:", n_population, "\n")
  cat("Number of PUMAs:", n_states, "\n")
  
  # Ensure age factor is properly set
  mi_data <- mi_data %>%
    mutate(
      age_category = factor(age_category,
                            levels = c("18-29", "30-44", "45-59", "60-74", "75+"))
    )
  
  # Get number of age categories
  age_categories <- levels(mi_data$age_category)
  n_age_categories <- length(age_categories)
  
  cat("\nAge categories:", n_age_categories, "\n")
  
  # Design matrix - using ONLY age
  X <- model.matrix(~ age_category, data = mi_data)
  cat("\nDesign matrix dimensions:", dim(X), "\n")
  cat("Column names:", colnames(X), "\n")
  
  # ========================================
  # STEP 1: GENERATE Y2 (Auxiliary Variable)
  # ========================================
  
  cat("\n--- Generating Y2 (auxiliary variable) ---\n")
  
  # Parameters for Y2 model - RANDOM DRAWS for age effects only
  gamma_intercept <- -1.0
  gamma_age <- c(0, rnorm(n_age_categories - 1, mean = 0, sd = 0.4))
  
  # Print drawn parameters
  cat("\nDrawn Y2 parameters:\n")
  cat("gamma_intercept:", round(gamma_intercept, 3), "\n")
  cat("gamma_age:", round(gamma_age, 3), "\n")
  
  # Combine gamma coefficients (intercept + age effects)
  gamma_all <- c(gamma_intercept, gamma_age[-1])
  names(gamma_all) <- colnames(X)
  
  # Random effects for Y2
  sigma_v <- 0.2
  v_j <- rnorm(n_states, 0, sigma_v)
  names(v_j) <- unique_states
  
  # Linear predictor for Y2
  eta_Y2 <- as.vector(X %*% gamma_all)
  
  # Add state random effects
  for (i in 1:n_population) {
    state <- mi_data$ST[i]
    eta_Y2[i] <- eta_Y2[i] + v_j[state]
  }
  
  # Generate Y2
  prob_Y2 <- inv_logit(eta_Y2)
  Y_2 <- rbinom(n_population, 1, prob_Y2)
  
  cat("\nY2 prevalence:", round(mean(Y_2), 3), "\n")
  
  # ========================================
  # STEP 2: GENERATE Y1 (Outcome depending on Y2)
  # ========================================
  
  cat("\n--- Generating Y1 (outcome) ---\n")
  
  # Parameters for Y1 model - RANDOM DRAWS for age effects only
  beta_intercept <- -1.5
  beta_age <- c(0, rnorm(n_age_categories - 1, mean = 0, sd = 0.5))
  
  # KEY PARAMETER: Effect of Y2 on Y1 (ensure positive)
  alpha <- abs(rnorm(1, mean = 1.5, sd = 0.3))
  
  # Print drawn parameters
  cat("\nDrawn Y1 parameters:\n")
  cat("beta_intercept:", round(beta_intercept, 3), "\n")
  cat("beta_age:", round(beta_age, 3), "\n")
  cat("alpha (Y2->Y1 effect):", round(alpha, 3), "\n")
  
  # Combine beta coefficients (intercept + age effects)
  beta_all <- c(beta_intercept, beta_age[-1])
  names(beta_all) <- colnames(X)
  
  # Random effects for Y1
  sigma_u <- 0.15
  u_j <- rnorm(n_states, 0, sigma_u)
  names(u_j) <- unique_states
  
  # Create extended design matrix including Y2
  X_extended <- cbind(X, Y2 = Y_2)
  beta_extended <- c(beta_all, alpha)
  
  # Linear predictor for Y1
  eta_Y1 <- as.vector(X_extended %*% beta_extended)
  
  # Add state random effects
  for (i in 1:n_population) {
    state <- mi_data$ST[i]
    eta_Y1[i] <- eta_Y1[i] + u_j[state]
  }
  
  # Generate Y1
  prob_Y1 <- inv_logit(eta_Y1)
  Y_1 <- rbinom(n_population, 1, prob_Y1)
  
  cat("Y1 prevalence:", round(mean(Y_1), 3), "\n")
  
  # ========================================
  # STEP 3: ANALYZE RELATIONSHIP
  # ========================================
  
  cat("\n--- Relationship Analysis ---\n")
  
  # Cross-tabulation
  cross_tab <- table(Y_1, Y_2)
  cat("\nCross-tabulation of Y1 and Y2:\n")
  print(cross_tab)
  
  # Calculate implied "sensitivity" and "specificity"
  # These are NOT measurement error parameters, but describe the association
  implied_sens <- sum(Y_1 == 1 & Y_2 == 1) / sum(Y_1 == 1)
  implied_spec <- sum(Y_1 == 0 & Y_2 == 0) / sum(Y_1 == 0)
  
  cat("\nImplied association measures:\n")
  cat("P(Y2=1|Y1=1):", round(implied_sens, 3), "\n")
  cat("P(Y2=0|Y1=0):", round(implied_spec, 3), "\n")
  cat("Correlation(Y1,Y2):", round(cor(Y_1, Y_2), 3), "\n")
  
  # By age (the only demographic)
  age_summary <- mi_data %>%
    mutate(Y_1 = Y_1, Y_2 = Y_2) %>%
    group_by(age_category) %>%
    summarise(
      n = n(),
      Y1_prev = mean(Y_1),
      Y2_prev = mean(Y_2),
      .groups = "drop"
    )
  cat("\nY1 and Y2 prevalence by age:\n")
  print(age_summary)
  
  # ========================================
  # STEP 4: PREPARE OUTPUT
  # ========================================
  
  # Add outcomes to data
  mi_data_with_outcomes <- mi_data %>%
    mutate(
      Y_1 = Y_1,
      Y_2 = Y_2
    )
  
  # Calculate true values by state
  true_values <- mi_data_with_outcomes %>%
    group_by(ST) %>%
    summarise(mean_Y1 = mean(Y_1), .groups = "drop") %>%
    deframe()
  
  # Add overall
  true_values <- c(true_values, overall = mean(Y_1))
  
  # Store all parameters
  true_params <- list(
    # Y2 model parameters
    gamma = gamma_all,
    gamma_intercept = gamma_intercept,
    gamma_age = gamma_age,
    v_j = v_j,
    sigma_v = sigma_v,
    
    # Y1 model parameters
    beta = beta_all,
    beta_intercept = beta_intercept,
    beta_age = beta_age,
    alpha = alpha,  # KEY: Y2 -> Y1 effect
    u_j = u_j,
    sigma_u = sigma_u,
    
    # For comparison with measurement error interpretation
    implied_sens = implied_sens,
    implied_spec = implied_spec,
    
    # Store category information
    n_age_categories = n_age_categories,
    age_categories = age_categories
  )
  
  cat("\n=== GENERATION COMPLETE ===\n")
  cat("Key parameter α (Y2->Y1 effect):", round(alpha, 3), "\n")
  cat("This is fundamentally different from measurement error!\n")
  cat("Y2 is a predictor/covariate, not a noisy measurement of Y1.\n")
  
  return(list(
    data = mi_data_with_outcomes,
    true_values = true_values,
    true_params = true_params,
    X = X
  ))
}

# ========================================
# USAGE EXAMPLE
# ========================================

# Load your Michigan data (assumes it's already loaded)
# mi_analysis_data <- readRDS("mi_pums_analysis_data.rds")

# Generate outcomes with age-only model
predictive_result_age <- generate_predictive_outcomes_age_only(mi_analysis_data, seed = 2025)

# Extract the data with outcomes
pop_data_predictive_age <- predictive_result_age$data
true_values_predictive_age <- predictive_result_age$true_values
true_params_predictive_age <- predictive_result_age$true_params

# Save for use with new Stan model
saveRDS(pop_data_predictive_age, "mi_population_predictive_age_only.rds")
saveRDS(true_values_predictive_age, "mi_true_values_predictive_age_only.rds")
saveRDS(true_params_predictive_age, "mi_true_params_predictive_age_only.rds")
