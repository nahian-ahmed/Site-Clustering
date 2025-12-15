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
plot_indiv <- FALSE          # Set to TRUE to generate individual plots per config
KEEP_BEST_FRACTION <- 1.0    # Filter to keep the top 100%

output_dir <- file.path("simulation_experiments", "output")
output_plot_dir <- file.path(output_dir, "plots")
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)

# Load data
true_params_df <- read.delim(file.path("config", "simulation_parameters.csv"), sep = ",", header = TRUE)
est_params_df  <- read.delim(file.path(output_dir ,"estimated_parameters.csv"), sep = ",", header = TRUE)
pred_perf_df   <- read.delim(file.path(output_dir ,"predictive_performance.csv"), sep = ",", header = TRUE)

# Grid Plot Configuration
grid_reference_methods <- c("0.125-kmSq", "clustGeo-50-80", "DBSC")
grid_param_sets <- c(1, 2, 3)

# Add a param_set index to true_params_df to merge easily
true_params_df$param_set <- 1:nrow(true_params_df)

# --- COLOR & FACTOR SETUP ---

# User specified order and colors
alg_all_o <- c(
    "reference",          
    "1to10",
    "2to10",
    "2to10-sameObs",
    "1-kmSq",
    "lat-long",
    "rounded-4",
    "SVS",
    "1-per-UL",           
    "DBSC",
    "BayesOptClustGeo"
)

colors <- c("red", "navy", "cyan", "pink", "green", "brown", "purple", "yellow", "blue", "darkgrey", "forestgreen")

# Create named color vector
color_map <- setNames(colors, alg_all_o)

# Identify methods to exclude from plotting
methods_to_exclude <- c("1-per-UL", "SVS")

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
  mutate(error = true_value - est_value) 

# -----------------------------------------------------------------------------
# 3. Filtering Simulations
# -----------------------------------------------------------------------------

cat(sprintf("Filtering simulations: Keeping top %.0f%% based on parameter closeness (MSE)...\n", KEEP_BEST_FRACTION * 100))

# A. Calculate Mean Squared Error (MSE) for each simulation run
sim_performance <- plot_data_full %>%
  group_by(reference_method, param_set, comparison_method, sim_num) %>%
  summarise(
    sim_mse = mean(error^2, na.rm = TRUE),
    valid_params = sum(!is.na(error)),
    total_params_count = n(),
    .groups = "drop"
  ) %>%
  filter(valid_params == total_params_count)

# B. Identify which simulations to keep
valid_sims_to_keep <- sim_performance %>%
  group_by(reference_method, param_set, comparison_method) %>%
  arrange(sim_mse) %>%
  mutate(
    rank = row_number(),
    total_runs = n(),
    cutoff_rank = ceiling(total_runs * KEEP_BEST_FRACTION)
  ) %>%
  filter(rank <= cutoff_rank) %>%
  select(reference_method, param_set, comparison_method, sim_num)

# C. Filter the main dataset
plot_data <- plot_data_full %>%
  inner_join(valid_sims_to_keep, by = c("reference_method", "param_set", "comparison_method", "sim_num"))

cat(sprintf("  - Original data points: %d\n", nrow(plot_data_full)))
cat(sprintf("  - Filtered data points: %d\n", nrow(plot_data)))

# Filter predictive performance data
pred_perf_df_filtered <- pred_perf_df %>%
  inner_join(valid_sims_to_keep, by = c("reference_method", "param_set", "comparison_method", "sim_num"))

# -----------------------------------------------------------------------------
# 4. Final Formatting
# -----------------------------------------------------------------------------

state_params <- c("state_intercept", "elevation", "TCB", "TCG", "TCW", "TCA")
obs_params   <- c("obs_intercept", "day_of_year", "time_observations_started", 
                  "duration_minutes", "effort_distance_km", "number_observers")

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
# 5. Generate 3x3 Summary Grid Plots
# -----------------------------------------------------------------------------

prepare_grid_data <- function(df, refs, params, exclude_methods, param_info_df, method_levels) {
  
  # 1. Base Filter
  base_df <- df %>%
    filter(reference_method %in% refs, 
           param_set %in% params,
           !comparison_method %in% exclude_methods)
  
  # 2. Add Beta_0 info for Column Labels
  base_df <- base_df %>%
    left_join(param_info_df %>% select(param_set, state_intercept), by = "param_set") %>%
    mutate(
      param_col_label = sprintf("Parameter~Set~%d~(beta[0] == %s)", 
                                param_set, 
                                state_intercept),
      ref_row_label = factor(paste("Reference =", reference_method),
                             levels = paste("Reference =", refs))
    )
  
  # 3. SMART DUPLICATION / RENAMING LOGIC
  df_refs <- base_df %>%
    filter(comparison_method == reference_method) %>%
    mutate(comparison_method = "reference")
  
  df_comps <- base_df %>%
    filter(comparison_method %in% method_levels, 
           comparison_method != "reference") 
  
  final_df <- bind_rows(df_refs, df_comps)
  
  # Set Factor Levels
  final_df$comparison_method <- factor(final_df$comparison_method, levels = method_levels)
  
  return(final_df)
}

# --- A. Parameter MAE Grid Plot (mae_grid.png) ---
cat("Generating 3x3 Parameter MAE Grid Plot (mae_grid.png)...\n")

mae_per_sim <- plot_data %>%
  group_by(reference_method, param_set, comparison_method, sim_num) %>%
  summarise(mae = mean(abs(error), na.rm = TRUE), .groups = "drop")

grid_mae_data <- prepare_grid_data(mae_per_sim, grid_reference_methods, grid_param_sets, methods_to_exclude, true_params_df, alg_all_o)

p_mae_grid <- ggplot(grid_mae_data, aes(x = comparison_method, y = mae, fill = comparison_method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  # Fixed scale ensures uniform y-axis across all subplots (removed scales="free_y")
  facet_grid(ref_row_label ~ param_col_label, labeller = labeller(.cols = label_parsed, .rows = label_value)) +
  scale_fill_manual(values = color_map) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(
    title = "Parameter Estimates: Mean Absolute Error (MAE)",
    # subtitle = "Rows = Reference Clustering | Columns = Parameter Sets",
    x = "Clustering Method",
    y = "Parameter MAE"
  )

ggsave(file.path(output_plot_dir, "mae_grid.png"), plot = p_mae_grid, width = 12, height = 10, dpi = 300)


# --- B. AUC Grid Plot (auc_grid.png) ---
cat("Generating 3x3 AUC Grid Plot (auc_grid.png)...\n")

grid_auc_data <- prepare_grid_data(pred_perf_df_filtered, grid_reference_methods, grid_param_sets, methods_to_exclude, true_params_df, alg_all_o)

p_auc_grid <- ggplot(grid_auc_data, aes(x = comparison_method, y = auc, fill = comparison_method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(ref_row_label ~ param_col_label, labeller = labeller(.cols = label_parsed, .rows = label_value)) +
  scale_fill_manual(values = color_map) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(
    title = "Predictive Performance: AUC",
    # subtitle = "Rows = Reference Clustering | Columns = Parameter Sets",
    x = "Clustering Method",
    y = "AUC"
  )

ggsave(file.path(output_plot_dir, "auc_grid.png"), plot = p_auc_grid, width = 12, height = 10, dpi = 300)


# --- C. AUPRC Grid Plot (auprc_grid.png) ---
cat("Generating 3x3 AUPRC Grid Plot (auprc_grid.png)...\n")

p_auprc_grid <- ggplot(grid_auc_data, aes(x = comparison_method, y = auprc, fill = comparison_method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(ref_row_label ~ param_col_label, labeller = labeller(.cols = label_parsed, .rows = label_value)) + 
  scale_fill_manual(values = color_map) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(
    title = "Predictive Performance: AUPRC",
    # subtitle = "Rows = Reference Clustering | Columns = Parameter Sets",
    x = "Clustering Method",
    y = "AUPRC"
  )

ggsave(file.path(output_plot_dir, "auprc_grid.png"), plot = p_auprc_grid, width = 12, height = 10, dpi = 300)


# -----------------------------------------------------------------------------
# 6. Aggregated Single Plots (Summarizing everything)
# -----------------------------------------------------------------------------

prepare_agg_data <- function(df, exclude_methods) {
  df %>%
    filter(!comparison_method %in% exclude_methods) %>%
    mutate(
      comparison_method = ifelse(comparison_method == reference_method, 
                                 "reference", 
                                 comparison_method)
    )
}

cat("Generating Aggregated Plots (mae.png, auc.png, auprc.png)...\n")

# --- MAE Aggregated (Sorted Descending) ---
agg_mae_data <- prepare_agg_data(mae_per_sim, methods_to_exclude)

mae_means <- agg_mae_data %>%
  group_by(comparison_method) %>%
  summarise(mean_val = mean(mae, na.rm = TRUE), .groups="drop") %>%
  arrange(desc(mean_val)) # Descending

agg_mae_data$comparison_method <- factor(agg_mae_data$comparison_method, levels = mae_means$comparison_method)

p_mae <- ggplot(agg_mae_data, aes(x = comparison_method, y = mae, fill = comparison_method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "darkred", color = "white") + # Add Means
  scale_fill_manual(values = color_map) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Aggregated Parameter MAE",
    # subtitle = "Across all Parameter Sets and Reference Clusterings (Sorted by Mean Descending)",
    x = "Clustering Method",
    y = "Parameter MAE"
  )
ggsave(file.path(output_plot_dir, "mae.png"), plot = p_mae, width = 6, height = 6, dpi = 300)


# --- AUC Aggregated (Sorted Ascending) ---
agg_auc_data <- prepare_agg_data(pred_perf_df_filtered, methods_to_exclude)

auc_means <- agg_auc_data %>%
  group_by(comparison_method) %>%
  summarise(mean_val = mean(auc, na.rm = TRUE), .groups="drop") %>%
  arrange(mean_val) # Ascending

agg_auc_data$comparison_method <- factor(agg_auc_data$comparison_method, levels = auc_means$comparison_method)

p_auc <- ggplot(agg_auc_data, aes(x = comparison_method, y = auc, fill = comparison_method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "darkred", color = "white") + # Add Means
  scale_fill_manual(values = color_map) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Aggregated AUC",
    # subtitle = "Across all Parameter Sets and Reference Clusterings (Sorted by Mean Ascending)",
    x = "Clustering Method",
    y = "AUC"
  )
ggsave(file.path(output_plot_dir, "auc.png"), plot = p_auc, width = 6, height = 6, dpi = 300)


# --- AUPRC Aggregated (Sorted Ascending) ---
agg_auprc_data <- prepare_agg_data(pred_perf_df_filtered, methods_to_exclude)

auprc_means <- agg_auprc_data %>%
  group_by(comparison_method) %>%
  summarise(mean_val = mean(auprc, na.rm = TRUE), .groups="drop") %>%
  arrange(mean_val) # Ascending

agg_auprc_data$comparison_method <- factor(agg_auprc_data$comparison_method, levels = auprc_means$comparison_method)

p_auprc <- ggplot(agg_auprc_data, aes(x = comparison_method, y = auprc, fill = comparison_method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "darkred", color = "white") + # Add Means
  scale_fill_manual(values = color_map) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Aggregated AUPRC",
    # subtitle = "Across all Parameter Sets and Reference Clusterings (Sorted by Mean Ascending)",
    x = "Clustering Method",
    y = "AUPRC"
  )
ggsave(file.path(output_plot_dir, "auprc.png"), plot = p_auprc, width = 6, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 7. Existing Individual Plots (Conditional)
# -----------------------------------------------------------------------------

if (plot_indiv) {
  cat("\nGenerating individual simulation plots (plot_indiv = TRUE)...\n")
  
  # --- Helper Functions ---
  make_error_plot <- function(ref_meth, p_set, data) {
    sub_data <- data %>% filter(reference_method == ref_meth, param_set == p_set)
    if(nrow(sub_data) == 0) return(NULL)
    
    p_state <- sub_data %>%
      filter(parameter %in% state_params) %>%
      ggplot(aes(x = comparison_method, y = error)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~parameter_label, nrow = 1, scales = "free_y") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = NULL, y = "Error (True - Est)")
    
    p_obs <- sub_data %>%
      filter(parameter %in% obs_params) %>%
      ggplot(aes(x = comparison_method, y = error)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~parameter_label, nrow = 1, scales = "free_y") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Clustering Method", y = "Error (True - Est)")
    
    return(p_state / p_obs)
  }

  make_pred_plot <- function(ref_meth, p_set, data, metric) {
    sub_data <- data %>% filter(reference_method == ref_meth, param_set == p_set)
    if(nrow(sub_data) == 0) return(NULL)
    y_lab <- toupper(metric)
    
    ggplot(sub_data, aes(x = comparison_method, y = .data[[metric]])) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.7, fill = "lightblue") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Clustering Method", y = y_lab,
           title = sprintf("Predictive Performance: %s (Ref: %s, Set: %d)", y_lab, ref_meth, p_set))
  }

  # --- Individual Plot Loop ---
  combinations <- unique(plot_data[, c("reference_method", "param_set")])
  plot_data_filtered_methods <- plot_data %>% filter(!comparison_method %in% methods_to_exclude)
  pred_perf_df_filtered_methods <- pred_perf_df_filtered %>% filter(!comparison_method %in% methods_to_exclude)

  for(i in 1:nrow(combinations)) {
    ref <- combinations$reference_method[i]
    par <- combinations$param_set[i]
    
    cat(sprintf("  Processing Ref=%s, ParamSet=%d...\n", ref, par))
    
    # Error Plots
    p <- make_error_plot(ref, par, plot_data)
    if(!is.null(p)) ggsave(file.path(output_plot_dir, sprintf("par_error_ref=%s_par=%d.png", ref, par)), plot = p, width = 14, height = 8, dpi = 300)
    
    p_filt <- make_error_plot(ref, par, plot_data_filtered_methods)
    if(!is.null(p_filt)) ggsave(file.path(output_plot_dir, sprintf("par_error_ref=%s_par=%d_no_SVS_UL.png", ref, par)), plot = p_filt, width = 14, height = 8, dpi = 300)
    
    # Prediction Plots
    p_auc <- make_pred_plot(ref, par, pred_perf_df, "auc")
    if(!is.null(p_auc)) ggsave(file.path(output_plot_dir, sprintf("auc_ref=%s_par=%d.png", ref, par)), plot = p_auc, width = 8, height = 6, dpi = 300)
    
    p_auprc <- make_pred_plot(ref, par, pred_perf_df, "auprc")
    if(!is.null(p_auprc)) ggsave(file.path(output_plot_dir, sprintf("auprc_ref=%s_par=%d.png", ref, par)), plot = p_auprc, width = 8, height = 6, dpi = 300)
  }
}

cat("All plots generated successfully.\n")