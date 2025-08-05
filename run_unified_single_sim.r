# ============================================
# UNIFIED SINGLE SIMULATION ANALYSIS
# For Models 1-4 with appropriate Stan files
# ============================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)
library(survey)

# Source required functions
source("generate_unified_dgp.R")  # The unified DGP for models 1-4

model_num = 4
overlap_pct = 0.11
seed = 17
alpha = 1.5
sigma_gamma = 2.5
demo_y2_interaction_sd = 1.0
N_ps = 1100; N_nps = 2500

# Function to run single simulation for any model
run_unified_single_sim <- function(model_num, 
                                   mi_analysis_data,
                                   seed = 17,
                                   alpha = 1.5,
                                   sigma_gamma = 2.5,
                                   demo_y2_interaction_sd = 1.0,
                                   N_ps = 1100,
                                   N_nps = 2500,
                                   overlap_pct = 0.11) {
  
  # ============================================
  # 1. GENERATE POPULATION
  # ============================================
  
  cat("\n=== GENERATING POPULATION FOR MODEL", model_num, "===\n")
  
  # Generate population using unified DGP
  dgp_result <- generate_unified_dgp(
    mi_data = mi_analysis_data,
    model = model_num,
    seed = seed,
    alpha = alpha,
    sigma_gamma = sigma_gamma,
    demo_y2_interaction_sd = demo_y2_interaction_sd
  )
  
  pop_data <- dgp_result$data
  true_values <- dgp_result$true_values
  true_params <- dgp_result$parameters
  
  # ============================================
  # 1b. CALCULATE TRUE TOTAL Y2 EFFECTS
  # ============================================
  
  cat("\n=== CALCULATING TRUE TOTAL Y2 EFFECTS ===\n")
  
  all_states <- sort(unique(pop_data$ST))
  true_total_y2_effects <- vector("numeric", length(all_states))
  names(true_total_y2_effects) <- all_states
  
  # Calculate true total Y2 effects
  for (j in seq_along(all_states)) {
    area <- all_states[j]
    area_data <- pop_data %>% filter(ST == area)
    X_area <- model.matrix(~ race_category + gender_binary + age_category, data = area_data)
    
    if (model_num == 1 || model_num == 2) {
      # Models 1 & 2: Y2 effect is just alpha
      true_total_y2_effects[area] <- alpha
      
    } else if (model_num == 3) {
      # Model 3: Y2 effect is alpha + X*delta (area-specific due to demographics)
      individual_effects <- alpha + as.vector(X_area %*% true_params$delta)
      true_total_y2_effects[area] <- mean(individual_effects)
      
    } else if (model_num == 4) {
      # Model 4: Y2 effect is alpha + gamma_j + X*delta
      individual_effects <- alpha + true_params$gamma_j[area] + 
        as.vector(X_area %*% true_params$delta)
      true_total_y2_effects[area] <- mean(individual_effects)
    }
  }
  
  cat("True total Y2 effects by area:\n")
  cat("Range: [", round(min(true_total_y2_effects), 3), ",", 
      round(max(true_total_y2_effects), 3), "]\n")
  
  # ============================================
  # 2. CREATE SAMPLES
  # ============================================
  
  cat("\n=== CREATING SAMPLES ===\n")
  
  samples <- create_overlap_samples_artificial(
    pop_data = pop_data,
    N_ps = N_ps,
    N_nps = N_nps,
    overlap_pct = overlap_pct,
    sampling_scheme = "SRS",
    use_mnar = FALSE,
    seed = seed
  )
  
  ps_sample <- samples$ps_data
  nps_sample <- samples$nps_data
  overlap_sample <- samples$overlap_data
  stan_data <- samples$stan_data
  
  # ============================================
  # 3. SELECT AND FIT MODELS
  # ============================================
  
  cat("\n=== SELECTING APPROPRIATE MODELS ===\n")
  
  # Select OG model based on model number
  og_stan_file <- switch(model_num,
                         "og_model_1.stan",  # Model 1
                         "og_model_2.stan",  # Model 2
                         "og_model_3.stan",  # Model 3
                         "og_enhanced_no_area_effects.stan"  # Model 4
  )
  
  # Select MRP model based on model number
  mrp_stan_file <- switch(model_num,
                          "MRP_predictive_no_area_effects.stan",   # Model 1: no area effects
                          "MRP_predictive.stan",                   # Model 2: with area effects
                          "MRP_predictive.stan",                   # Model 3: with area effects
                          "MRP_predictive.stan"                    # Model 4: with area effects
                          # "MRP_predictive.stan"
  )
  
  cat("OG Model:", og_stan_file, "\n")
  cat("MRP Model:", mrp_stan_file, "\n")
  
  # Compile models
  og_model <- cmdstan_model(og_stan_file)
  mrp_model <- cmdstan_model(mrp_stan_file)
  
  # Fit OG
  cat("\nFitting OG model...\n")
  fit_og <- og_model$sample(
    data = stan_data,
    seed = 23,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 100
  )
  
  # Fit MRP
  cat("\nFitting MRP model...\n")
  fit_mrp <- mrp_model$sample(
    data = stan_data,
    seed = 31,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 100
  )
  
  # ============================================
  # 4. EXTRACT ESTIMATES
  # ============================================
  
  cat("\n=== EXTRACTING ESTIMATES ===\n")
  
  # Extract posterior draws
  posterior_probs <- fit_og$draws("prob_Y1_given_Y2", format = "df") |> as_tibble()
  posterior_nps <- fit_og$draws("group_mrp_estimate", format = "df") |> as_tibble()
  overall_nps <- fit_og$draws("overall_mrp_estimate", format = "draws_matrix") |> as.numeric()
  
  # Extract estimated total Y2 effects from OG model (if available)
  if (model_num %in% c(3, 4)) {
    avg_total_y2_draws <- fit_og$draws("avg_total_y2_effect_by_area", format = "draws_matrix")
    estimated_total_y2_effects <- apply(avg_total_y2_draws, 2, mean)
    y2_effect_lower <- apply(avg_total_y2_draws, 2, quantile, 0.025)
    y2_effect_upper <- apply(avg_total_y2_draws, 2, quantile, 0.975)
  } else {
    # For models 1 and 2, Y2 effect is constant
    estimated_total_y2_effects <- rep(alpha, length(all_states))
    y2_effect_lower <- rep(0, length(all_states))
    y2_effect_upper <- rep(0, length(all_states))
  }
  
  posterior_mrp <- fit_mrp$draws("group_mrp_estimate", format = "df") |> as_tibble()
  overall_mrp <- fit_mrp$draws("overall_mrp_estimate", format = "draws_matrix") |> as.numeric()
  
  num_draws <- nrow(posterior_probs)
  
  # ============================================
  # 5. CALCULATE PS-BASED ESTIMATES
  # ============================================
  
  cat("\n=== CALCULATING PS-BASED ESTIMATES ===\n")
  
  group_ps_estimates <- lapply(seq_len(num_draws), function(m) {
    Y_probs <- posterior_probs[m, 1:nrow(ps_sample)] |> unlist()
    ps_with_probs <- ps_sample |> mutate(Y_val = Y_probs)
    overlap_with_Y1 <- overlap_sample |> mutate(Y_val = as.numeric(Y_1))
    combined_data <- bind_rows(ps_with_probs, overlap_with_Y1)
    dsgn <- create_survey_design(combined_data, sampling_scheme = "SRS")
    svyby(~Y_val, ~ST, dsgn, svymean, keep.var = TRUE, deff = FALSE)
  })
  
  # Extract point estimates
  ps_estimates <- sapply(all_states, function(st) {
    vals <- sapply(group_ps_estimates, function(x) x[x$ST == st, "Y_val"])
    mean(vals)
  })
  
  nps_estimates <- sapply(seq_along(all_states), function(j) {
    nps_cols <- posterior_nps |> select(starts_with("group_mrp_estimate"))
    mean(nps_cols[[j]])
  })
  
  mrp_estimates <- sapply(seq_along(all_states), function(j) {
    mrp_cols <- posterior_mrp |> select(starts_with("group_mrp_estimate"))
    mean(mrp_cols[[j]])
  })
  
  # ============================================
  # 6. CREATE COMPARISON DATA
  # ============================================
  
  # Get area analysis from true_params
  area_analysis <- true_params$area_analysis
  
  comparison_df <- tibble(
    area = all_states,
    true_value = unname(true_values[all_states]),
    og_ps = ps_estimates,
    og_nps = nps_estimates,
    mrp = mrp_estimates,
    true_total_y2_effect = true_total_y2_effects,
    est_total_y2_effect = estimated_total_y2_effects,
    est_y2_lower = y2_effect_lower,
    est_y2_upper = y2_effect_upper
  ) |>
    mutate(
      bias_ps = og_ps - true_value,
      bias_nps = og_nps - true_value,
      bias_mrp = mrp - true_value,
      y2_effect_error = est_total_y2_effect - true_total_y2_effect
    ) |>
    left_join(
      area_analysis |> select(ST, Y2_prev),
      by = c("area" = "ST")
    ) |>
    mutate(
      y2_impact = Y2_prev * true_total_y2_effect
    )
  
  
  # ============================================
  # 7. POSTERIOR PREDICTIVE CHECKS
  # ============================================
  
  cat("\n=== GENERATING PPCs ===\n")
  
  # Get true Y1 for PS sample
  y1_data <- left_join(
    ps_sample |> select(person_id, ST, age_category, race_category, gender_binary),
    pop_data |> select(person_id, Y_1),
    by = "person_id"
  )
  y1_true_ps <- y1_data$Y_1
  
  # OG PPCs
  y1_rep_og <- fit_og$draws("Y1_ps_imputed", format = "matrix")
  
  # MRP PPCs - Uses NPS sample 
  y1_rep_mrp <- fit_mrp$draws("Y1_nps_rep", format = "matrix")
  y1_true_nps <- nps_sample$Y_1
  
  # ============================================
  # 8. CREATE PLOTS
  # ============================================
  
  cat("\n=== CREATING PLOTS ===\n")
  
  # 1. Point estimates vs truth
  p_estimates <- comparison_df |>
    pivot_longer(cols = c(og_ps, og_nps, mrp), 
                 names_to = "method", 
                 values_to = "estimate") |>
    mutate(method = case_when(
      method == "og_ps" ~ "OG-PS",
      method == "og_nps" ~ "OG-NPS",
      method == "mrp" ~ "MRP"
    )) |>
    ggplot(aes(x = true_value, y = estimate, color = method)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_manual(values = c("OG-PS" = "#E41A1C", 
                                  "OG-NPS" = "#377EB8", 
                                  "MRP" = "#984EA3")) +
    labs(title = "Area-level Estimates vs Truth",
         x = "True Prevalence", 
         y = "Estimated Prevalence") +
    theme_minimal() +
    facet_wrap(~method, ncol = 3)
  
  # 2. Bias patterns by COMPLETE Y2 impact
  p_bias <- comparison_df |>
    pivot_longer(cols = c(bias_nps, bias_mrp), 
                 names_to = "method", 
                 values_to = "bias") |>
    mutate(method = ifelse(method == "bias_nps", "OG-NPS", "MRP")) |>
    ggplot(aes(x = y2_impact, y = bias, color = method)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("OG-NPS" = "#377EB8", "MRP" = "#984EA3")) +
    labs(title = "Bias Pattern by Complete Y2 Impact",
         subtitle = "Y2 Impact = Y2 Prevalence × Total Y2 Effect (α + γ_j + avg(X*δ))",
         x = "Y2 Impact", 
         y = "Bias") +
    theme_minimal()
  
  # 3. PPCs - Overall
  ppc_overall <- (
    bayesplot::ppc_bars(y1_true_ps, y1_rep_og) +
      labs(title = "PPC Overall: OG", subtitle = "True vs Predicted Y1")
  ) + (
    bayesplot::ppc_bars(y1_true_nps, y1_rep_mrp) +
      labs(title = "PPC Overall: MRP", subtitle = "True vs Predicted Y1")
  )
  
  # 4. PPCs - By demographics
  ppc_age <- (
    bayesplot::ppc_bars_grouped(y1_true_ps, y1_rep_og, 
                                group = ps_sample$age_category) +
      labs(title = "PPC by Age: OG")
  ) + (
    bayesplot::ppc_bars_grouped(y1_true_nps, y1_rep_mrp, 
                                group = nps_sample$age_category) +
      labs(title = "PPC by Age: MRP")
  )
  
  ppc_area <- (
    bayesplot::ppc_bars_grouped(y1_true_ps, y1_rep_og, 
                                group = ps_sample$ST) +
      labs(title = "PPC by Area: OG")
  ) | (
    bayesplot::ppc_bars_grouped(y1_true_nps, y1_rep_mrp, 
                                group = nps_sample$ST) +
      labs(title = "PPC by Area: MRP")
  )
  
  ppc_race <- (
    bayesplot::ppc_bars_grouped(y1_true_ps, y1_rep_og, 
                                group = ps_sample$race_category) +
      labs(title = "PPC by Race: OG")
  ) + (
    bayesplot::ppc_bars_grouped(y1_true_nps, y1_rep_mrp, 
                                group = nps_sample$race_category) +
      labs(title = "PPC by Race: MRP")
  )
  
  ppc_gender <- (
    bayesplot::ppc_bars_grouped(y1_true_ps, y1_rep_og, 
                                group = ps_sample$gender_binary) +
      labs(title = "PPC by Gender: OG")
  ) + (
    bayesplot::ppc_bars_grouped(y1_true_nps, y1_rep_mrp, 
                                group = nps_sample$gender_binary) +
      labs(title = "PPC by Gender: MRP")
  )
  
  # 5. Y2 effect comparison: True vs Estimated
  p_y2_recovery <- comparison_df |>
    ggplot(aes(x = true_total_y2_effect, y = est_total_y2_effect)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = "Y2 Effect Recovery: OG Model",
         subtitle = "Total Y2 Effect = α + γ_j + avg(X*δ)",
         x = "True Total Y2 Effect", 
         y = "Estimated Total Y2 Effect") +
    theme_minimal()
  
  # 6. Y2 effect visualization - showing ALL components with CIs
  p_y2_effects <- comparison_df |>
    arrange(true_total_y2_effect) |>
    mutate(area = factor(area, levels = area)) |>
    ggplot(aes(x = area)) +
    geom_col(aes(y = true_total_y2_effect, fill = "True Total"), 
             alpha = 0.7, position = "identity") +
    geom_point(aes(y = est_total_y2_effect, color = "Estimated"), 
               size = 3) +
    geom_errorbar(aes(ymin = est_y2_lower, ymax = est_y2_upper), 
                  width = 0.3, color = "darkred", alpha = 0.8) +
    geom_hline(yintercept = true_params$alpha, 
               linetype = "dashed", color = "red") +
    coord_flip() +
    scale_fill_manual(values = c("True Total" = "steelblue")) +
    scale_color_manual(values = c("Estimated" = "darkred")) +
    labs(title = "Total Y2 Effect by Area: True vs Estimated with 95% CIs",
         subtitle = paste("Red line = global effect α =", true_params$alpha, 
                          "\nBars: true effects, Points: estimates, Error bars: 95% credible intervals"),
         x = "Area", 
         y = "Total Y2 Effect") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # 7. Bias comparison
  bias_summary <- comparison_df |>
    summarise(
      `OG-PS` = sqrt(mean(bias_ps^2)),
      `OG-NPS` = sqrt(mean(bias_nps^2)),
      MRP = sqrt(mean(bias_mrp^2))
    ) |>
    pivot_longer(everything(), names_to = "Method", values_to = "RMSE")
  
  p_rmse <- bias_summary |>
    ggplot(aes(x = Method, y = RMSE, fill = Method)) +
    geom_col() +
    scale_fill_manual(values = c("OG-PS" = "#E41A1C", 
                                 "OG-NPS" = "#377EB8", 
                                 "MRP" = "#984EA3")) +
    labs(title = "Root Mean Square Error",
         y = "RMSE") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 8. Parameter recovery - expanded
  if (model_num %in% c(3, 4)) {
    param_recovery <- bayesplot::mcmc_recover_intervals(
      fit_og$draws(c("beta", "alpha",
                     "delta_raw",
                     "sigma_gamma")
      ),
      true = c(true_params$beta, true_params$alpha,
               true_params$delta[-1],
               true_params$sigma_gamma)
    ) + 
      labs(title = "Parameter Recovery: OG Model",
           subtitle = "True values (dots) vs 95% posterior intervals") +
      theme_minimal()
  } else if (model_num == 2) {
    param_recovery <- bayesplot::mcmc_recover_intervals(
      fit_og$draws(c("beta", "alpha",
                     "sigma_gamma")
      ),
      true = c(true_params$beta, true_params$alpha,
               true_params$sigma_gamma)
    ) + 
      labs(title = "Parameter Recovery: OG Model",
           subtitle = "True values (dots) vs 95% posterior intervals") +
      theme_minimal()
  } else if (model_num == 1) {
    param_recovery <- bayesplot::mcmc_recover_intervals(
      fit_og$draws(c("beta", "alpha")
      ),
      true = c(true_params$beta, true_params$alpha)
    ) + 
      labs(title = "Parameter Recovery: OG Model",
           subtitle = "True values (dots) vs 95% posterior intervals") +
      theme_minimal()
  }
  
  # ============================================
  # 9. COMBINE AND DISPLAY
  # ============================================
  
  # Main comparison
  # main_plot <- (p_estimates / p_bias) / (p_y2_recovery | p_rmse) +
  #   plot_annotation(
  #     title = "Single Simulation: OG vs MRP with Strong Y2 Interactions",
  #     # subtitle = "Complete Y2 effects",
  #     theme = theme(plot.title = element_text(size = 18, face = "bold"))
  #   )
  
  main_plot <- (p_estimates | p_rmse) +
    plot_annotation(
      title = paste0("Model ", model_num, " Single Simulation: OG vs MRP"),
      subtitle = true_params$model_desc
    ) & 
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  
  # PPC comparison
  ppc_plot1 <- ppc_overall / ppc_area +
    plot_layout(heights = c(2, 5)) +  # First plot is 1 unit, second is 2 units
    plot_annotation(
      title = "Posterior Predictive Checks: OG vs MRP",
      subtitle = "How well do models capture demographic patterns?",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  ppc_plot2 <- ppc_age / ppc_race / ppc_gender +
    plot_annotation(
      title = "Posterior Predictive Checks: OG vs MRP",
      subtitle = "How well do models capture demographic patterns?",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Y2 effect plot
  y2_plot <- p_y2_effects +
    plot_annotation(
      title = "Y2 Effect Structure",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # # Display plots
  # print(main_plot)
  # print(ppc_plot1)
  # print(ppc_plot2)
  # print(y2_plot)
  # print(param_recovery)
  
  # ============================================
  # 10. SUMMARY STATISTICS
  # ============================================
  
  cat("\n=== SUMMARY FOR SINGLE SIMULATION ===\n")
  
  cat("\nRMSE by Method:\n")
  print(bias_summary)
  
  cat("\nWorst 5 areas for MRP:\n")
  worst_mrp <- comparison_df |>
    arrange(desc(abs(bias_mrp))) |>
    slice_head(n = 5) |>
    select(area, true_value, mrp, bias_mrp, Y2_prev, true_total_y2_effect)
  print(worst_mrp)
  
  cat("\nBias correlation with COMPLETE Y2 impact:\n")
  cor_og <- with(comparison_df, cor(bias_nps, y2_impact))
  cor_mrp <- with(comparison_df, cor(bias_mrp, y2_impact))
  cat("OG-NPS:", round(cor_og, 3), "\n")
  cat("MRP:", round(cor_mrp, 3), "\n")
  
  cat("\nY2 Effect Recovery (OG):\n")
  y2_recovery_cor <- cor(comparison_df$true_total_y2_effect, 
                         comparison_df$est_total_y2_effect)
  y2_recovery_rmse <- sqrt(mean((comparison_df$true_total_y2_effect - 
                                   comparison_df$est_total_y2_effect)^2))
  cat("Correlation:", round(y2_recovery_cor, 3), "\n")
  cat("RMSE:", round(y2_recovery_rmse, 3), "\n")
  
  
  # ============================================
  # 11. SAVE RESULTS
  # ============================================
  
  output_prefix <- paste0("model", model_num)
  
  # Save plots
  ggsave(paste0(output_prefix, "_main_plot.png"), main_plot, 
         width = 14, height = 10, dpi = 150)
  ggsave(paste0(output_prefix, "_ppc_plot1.png"), ppc_plot1, 
         width = 14, height = 10, dpi = 150)
  ggsave(paste0(output_prefix, "_ppc_plot2.png"), ppc_plot2, 
         width = 14, height = 10, dpi = 150)
  ggsave(paste0(output_prefix, "_y2_plot.png"), y2_plot, 
         width = 14, height = 10, dpi = 150)
  ggsave(paste0(output_prefix, "_param_recovery.png"), param_recovery, 
         width = 14, height = 10, dpi = 150)
  
  
  
  
  # Save comparison data
  saveRDS(comparison_df, paste0(output_prefix, "_comparison_df.rds"))
  saveRDS(bias_summary, paste0(output_prefix, "_bias_summary.rds"))
  
  # Return results
  return(
    list(
      comparison_df = comparison_df,
      bias_summary = bias_summary,
      fits = list(og = fit_og, mrp = fit_mrp),
      true_params = true_params,
      plots = list(
        main = main_plot,
        ppc1 = ppc_plot1,
        ppc2 = ppc_plot2,
        y2 = y2_plot,
        param_rec = param_recovery
      )
    )
  )

}


# Save fits
saveRDS(list(og = fit_og, mrp = fit_mrp), paste0(output_prefix, "_fits.rds"))

# ============================================
# 9. PRINT SUMMARY
# ============================================

cat("\n=== SUMMARY FOR MODEL", model_num, "===\n")
cat(true_params$model_desc, "\n\n")

cat("RMSE by Method:\n")
print(bias_summary)

cat("\nWorst 5 areas for MRP:\n")
worst_mrp <- comparison_df |>
  arrange(desc(abs(bias_mrp))) |>
  slice_head(n = 5) |>
  select(area, true_value, mrp, bias_mrp, Y2_prev, true_total_y2_effect)
print(worst_mrp)

cat("\nBias correlation with Y2 impact:\n")
cor_og <- with(comparison_df, cor(bias_nps, y2_impact))
cor_mrp <- with(comparison_df, cor(bias_mrp, y2_impact))
cat("OG-NPS:", round(cor_og, 3), "\n")
cat("MRP:", round(cor_mrp, 3), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

# Return results
return(list(
  comparison_df = comparison_df,
  fits = list(og = fit_og, mrp = fit_mrp),
  true_params = true_params,
  plots = list(main = main_plot)
))

# ============================================
# RUN ALL MODELS
# ============================================

# Example usage:
if (TRUE) {  # Set to TRUE to run
  
  # Run all four models
  results <- list()
  
  for (model_num in 1:4) {
    cat("\n", rep("=", 80), "\n")
    cat("RUNNING MODEL", model_num, "\n")
    cat(rep("=", 80), "\n")
    
    results[[paste0("model", model_num)]] <- run_unified_single_sim(
      model_num = model_num,
      mi_analysis_data = mi_analysis_data,
      seed = 17,
      alpha = 1.5,
      sigma_gamma = 2.5,
      demo_y2_interaction_sd = 1.0
    )
  }
  
  # Save all results
  saveRDS(results, "all_models_single_sim_results.rds")
}
