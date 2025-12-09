# plot_simulations.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readr)
library(stringr)

# -----------------------------------------------------------------------------
# 1. Setup and Data Loading
# -----------------------------------------------------------------------------

# --- CONFIGURATION ---
KEEP_BEST_FRACTION <- 1.0 # Filter to keep the top 90% of simulations closest to truth

output_dir <- file.path("simulation_experiments", "output")
# Load data
true_params_df <- read.delim(file.path("config", "simulation_parameters.csv"), sep = ",", header = T)
est_params_df  <- read.delim(file.path(output_dir ,"estimated_parameters.csv"), sep = ",", header = T)
pred_perf_df   <- read.delim(file.path(output_dir ,"predictive_performance.csv"), sep = ",", header = T)


output_plot_dir <- file.path(output_dir, "plots")
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)

# Methods to exclude for the filtered plots
methods_to_exclude <- c("1-per-UL", "SVS")

# Add a param_set index to true_params_df to merge easily
true_params_df$param_set <- 1:nrow(true_params_df)

# -----------------------------------------------------------------------------
# 2. Data Processing for Parameter Errors
# -----------------------------------------------------------------------------

# Pivot true parameters to long format
true_long <- true_params_df %>%
  pivot_longer(
    cols = -param_set,
    names_to = "parameter",
    values_to = "true_value"
  )

# Pivot estimated parameters to long format
# Identifying parameter columns (excluding metadata columns)
meta_cols <- c("reference_method", "param_set", "sim_num", "comparison_method", "nll", "convergence")
param_cols <- setdiff(names(est_params_df), meta_cols)

est_long <- est_params_df %>%
  pivot_longer(
    cols = all_of(param_cols),
    names_to = "parameter",
    values_to = "est_value"
  )

# Merge true and estimated data
plot_data_full <- est_long %>%
  left_join(true_long, by = c("param_set", "parameter")) %>%
  mutate(error = true_value - est_value) # Error = True - Estimated

# -----------------------------------------------------------------------------
# 3. Filtering Top N% Simulations (Closest to Truth)
# -----------------------------------------------------------------------------

cat(sprintf("Filtering simulations: Keeping top %.0f%% based on parameter closeness (MSE)...\n", KEEP_BEST_FRACTION * 100))

# A. Calculate Mean Squared Error (MSE) for each simulation run
#    (Aggregating error across all parameters for that run)
sim_performance <- plot_data_full %>%
  group_by(reference_method, param_set, comparison_method, sim_num) %>%
  summarise(
    # Scalar metric for "closeness": Mean Squared Error across all params
    sim_mse = mean(error^2, na.rm = TRUE),
    # Count valid parameters to ensure we don't rank failed/NA runs as "good"
    valid_params = sum(!is.na(error)),
    total_params_count = n(),
    .groups = "drop"
  ) %>%
  # Exclude runs that have NAs (failed convergence)
  filter(valid_params == total_params_count)

# B. Identify which simulations to keep
valid_sims_to_keep <- sim_performance %>%
  group_by(reference_method, param_set, comparison_method) %>%
  arrange(sim_mse) %>%
  mutate(
    rank = row_number(),
    total_runs = n(),
    # Calculate cutoff rank (e.g., 90% of 10 runs = 9)
    cutoff_rank = ceiling(total_runs * KEEP_BEST_FRACTION)
  ) %>%
  filter(rank <= cutoff_rank) %>%
  select(reference_method, param_set, comparison_method, sim_num)

# C. Filter the main dataset
plot_data <- plot_data_full %>%
  inner_join(valid_sims_to_keep, by = c("reference_method", "param_set", "comparison_method", "sim_num"))

cat(sprintf("  - Original data points: %d\n", nrow(plot_data_full)))
cat(sprintf("  - Filtered data points: %d\n", nrow(plot_data)))


# -----------------------------------------------------------------------------
# 4. Final Formatting for Plotting
# -----------------------------------------------------------------------------

# Define parameter groups for the two rows
state_params <- c("state_intercept", "elevation", "TCB", "TCG", "TCW", "TCA")
obs_params   <- c("obs_intercept", "day_of_year", "time_observations_started", 
                  "duration_minutes", "effort_distance_km", "number_observers")

# Clean labels for plotting
label_map <- c(
  "state_intercept" = "State Intercept",
  "elevation" = "Elevation",
  "TCB" = "TCB", "TCG" = "TCG", "TCW" = "TCW", "TCA" = "TCA",
  "obs_intercept" = "Obs. Intercept",
  "day_of_year" = "Day of Year",
  "time_observations_started" = "Time Obs. Start",
  "duration_minutes" = "Duration (min)",
  "effort_distance_km" = "Effort Dist. (km)",
  "number_observers" = "Num. Observers"
)

plot_data$parameter_label <- factor(plot_data$parameter, 
                                    levels = c(state_params, obs_params),
                                    labels = label_map[c(state_params, obs_params)])

# -----------------------------------------------------------------------------
# 5. Function to Generate Parameter Error Plots
# -----------------------------------------------------------------------------

make_error_plot <- function(ref_meth, p_set, data) {
  
  # Filter data for current combination
  sub_data <- data %>%
    filter(reference_method == ref_meth, param_set == p_set)
  
  if(nrow(sub_data) == 0) return(NULL)
  
  # Prepare State Parameter Plots
  p_state <- sub_data %>%
    filter(parameter %in% state_params) %>%
    ggplot(aes(x = comparison_method, y = error)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~parameter_label, nrow = 1, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Error (True - Est)")
  
  # Prepare Observation Parameter Plots
  p_obs <- sub_data %>%
    filter(parameter %in% obs_params) %>%
    ggplot(aes(x = comparison_method, y = error)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~parameter_label, nrow = 1, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Clustering Method", y = "Error (True - Est)")
  
  # Combine rows using patchwork
  combined_plot <- p_state / p_obs
  
  return(combined_plot)
}

# -----------------------------------------------------------------------------
# 6. Loop to Create and Save Parameter Error Plots
# -----------------------------------------------------------------------------

combinations <- unique(plot_data[, c("reference_method", "param_set")])

for(i in 1:nrow(combinations)) {
  ref <- combinations$reference_method[i]
  par <- combinations$param_set[i]
  
  cat(sprintf("Generating error plots for Ref=%s, ParamSet=%d...\n", ref, par))
  
  # --- 1. FULL PLOT (Original) ---
  p <- make_error_plot(ref, par, plot_data)
  
  if(!is.null(p)) {
    filename <- sprintf("par_error_ref=%s_par=%d.png", ref, par)
    ggsave(file.path(output_plot_dir, filename), plot = p, width = 14, height = 8, dpi = 300)
  }

  # --- 2. FILTERED PLOT (No 1-per-UL or SVS) ---
  plot_data_filtered_methods <- plot_data %>% 
    filter(!comparison_method %in% methods_to_exclude)
    
  p_filtered <- make_error_plot(ref, par, plot_data_filtered_methods)
  
  if(!is.null(p_filtered)) {
    filename_filtered <- sprintf("par_error_ref=%s_par=%d_no_SVS_UL.png", ref, par)
    ggsave(file.path(output_plot_dir, filename_filtered), plot = p_filtered, width = 14, height = 8, dpi = 300)
  }
}

# -----------------------------------------------------------------------------
# 7. Function to Generate Predictive Performance Plots
# -----------------------------------------------------------------------------

make_pred_plot <- function(ref_meth, p_set, data, metric) {
  
  sub_data <- data %>%
    filter(reference_method == ref_meth, param_set == p_set)
  
  if(nrow(sub_data) == 0) return(NULL)
  
  # Determine y-label
  y_lab <- toupper(metric)
  
  p <- ggplot(sub_data, aes(x = comparison_method, y = .data[[metric]])) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7, fill = "lightblue") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Clustering Method", 
         y = y_lab,
         title = sprintf("Predictive Performance: %s (Ref: %s, Set: %d)", y_lab, ref, p_set))
  
  return(p)
}

# -----------------------------------------------------------------------------
# 8. Loop to Create and Save Predictive Performance Plots
# -----------------------------------------------------------------------------

# Use the same combinations from the predictive performance dataframe
pred_combinations <- unique(pred_perf_df[, c("reference_method", "param_set")])

for(i in 1:nrow(pred_combinations)) {
  ref <- pred_combinations$reference_method[i]
  par <- pred_combinations$param_set[i]
  
  cat(sprintf("Generating pred plots for Ref=%s, ParamSet=%d...\n", ref, par))
  
  # Create filtered dataframe
  pred_perf_df_filtered <- pred_perf_df %>% 
    filter(!comparison_method %in% methods_to_exclude)

  # --- AUC Plots ---
  # 1. Full
  p_auc <- make_pred_plot(ref, par, pred_perf_df, "auc")
  if(!is.null(p_auc)) {
    filename_auc <- sprintf("auc_ref=%s_par=%d.png", ref, par)
    ggsave(file.path(output_plot_dir, filename_auc), plot = p_auc, width = 8, height = 6, dpi = 300)
  }
  
  # 2. Filtered
  p_auc_filtered <- make_pred_plot(ref, par, pred_perf_df_filtered, "auc")
  if(!is.null(p_auc_filtered)) {
    filename_auc_filtered <- sprintf("auc_ref=%s_par=%d_no_SVS_UL.png", ref, par)
    ggsave(file.path(output_plot_dir, filename_auc_filtered), plot = p_auc_filtered, width = 8, height = 6, dpi = 300)
  }
  
  # --- AUPRC Plots ---
  # 1. Full
  p_auprc <- make_pred_plot(ref, par, pred_perf_df, "auprc")
  if(!is.null(p_auprc)) {
    filename_auprc <- sprintf("auprc_ref=%s_par=%d.png", ref, par)
    ggsave(file.path(output_plot_dir, filename_auprc), plot = p_auprc, width = 8, height = 6, dpi = 300)
  }
  
  # 2. Filtered
  p_auprc_filtered <- make_pred_plot(ref, par, pred_perf_df_filtered, "auprc")
  if(!is.null(p_auprc_filtered)) {
    filename_auprc_filtered <- sprintf("auprc_ref=%s_par=%d_no_SVS_UL.png", ref, par)
    ggsave(file.path(output_plot_dir, filename_auprc_filtered), plot = p_auprc_filtered, width = 8, height = 6, dpi = 300)
  }
}

cat("All plots generated successfully.\n")
