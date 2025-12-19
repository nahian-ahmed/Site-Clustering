# plot_simulation_ref.R

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. SETUP & DATA LOADING
# -----------------------------------------------------------------------------
output_dir <- file.path("simulation_experiments", "output")
plot_output_dir <- file.path(output_dir, "plots_ref")
if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)

# Load data and configurations
true_params_df <- read.csv(file.path("config", "simulation_parameters.csv")) #
sim_methods_df <- read.csv(file.path("config", "simulation_clusterings.csv")) #
est_params_df  <- read.csv(file.path(output_dir, "estimated_parameters.csv"))
pred_perf_df   <- read.csv(file.path(output_dir, "predictive_performance.csv"))

true_params_df$param_set <- seq_len(nrow(true_params_df))

# Set method order from simulation_clusterings.csv
method_order <- sim_methods_df$method

# Define parameter groups and order from simulation_parameters.csv headers
state_params <- c("state_intercept", "elevation", "TCB", "TCG", "TCW", "TCA")
obs_params   <- c("obs_intercept", "day_of_year", "time_observations_started", 
                  "duration_minutes", "effort_distance_km", "number_observers")

label_map <- c(
  "state_intercept" = "State Intercept", "elevation" = "Elevation",
  "TCB" = "TCB", "TCG" = "TCG", "TCW" = "TCW", "TCA" = "TCA",
  "obs_intercept" = "Obs. Intercept", "day_of_year" = "Day of Year",
  "time_observations_started" = "Time Start", "duration_minutes" = "Duration",
  "effort_distance_km" = "Effort Dist", "number_observers" = "Observers"
)

# -----------------------------------------------------------------------------
# 2. DATA PRE-PROCESSING
# -----------------------------------------------------------------------------
true_long <- true_params_df %>% 
  pivot_longer(-param_set, names_to = "parameter", values_to = "true_val")

est_long  <- est_params_df %>% 
  pivot_longer(cols = all_of(c(state_params, obs_params)), names_to = "parameter", values_to = "est_val") %>%
  left_join(true_long, by = c("param_set", "parameter")) %>%
  mutate(error = true_val - est_val)

# Enforce factor levels for parameters and methods
est_long$parameter <- factor(est_long$parameter, levels = c(state_params, obs_params))
est_long$comparison_method <- factor(est_long$comparison_method, levels = method_order)

# Prepare labels with Greek letters for TOP ROW only
facet_labeller <- true_params_df %>%
  select(param_set, state_intercept) %>%
  mutate(col_title = sprintf("Parameter~Set~%d~(beta[0] == %s)", param_set, state_intercept))

# -----------------------------------------------------------------------------
# 3. INDIVIDUAL PARAMETER ERROR PLOTS (TRUE - ESTIMATE)
# -----------------------------------------------------------------------------
combos <- unique(est_long[, c("reference_method", "param_set")])

for (i in 1:nrow(combos)) {
  ref <- combos$reference_method[i]
  ps  <- combos$param_set[i]
  
  sub_data <- est_long %>% filter(reference_method == ref, param_set == ps)
  
  p_state <- ggplot(sub_data %>% filter(parameter %in% state_params), 
                    aes(x = comparison_method, y = error, fill = comparison_method)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~parameter, scales = "free_y", labeller = as_labeller(label_map), nrow = 1) +
    theme_bw() + 
    theme(axis.text.x = element_blank(), legend.position = "none") +
    labs(y = "Error (True-Est)", x = NULL, title = sprintf("Ref: %s | Param Set: %d", ref, ps))

  p_obs <- ggplot(sub_data %>% filter(parameter %in% obs_params), 
                  aes(x = comparison_method, y = error, fill = comparison_method)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~parameter, scales = "free_y", labeller = as_labeller(label_map), nrow = 1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(y = "Error (True-Est)", x = "Method")

  ggsave(file.path(plot_output_dir, sprintf("ref=%s_par=%d.png", ref, ps)), p_state / p_obs, width = 12, height = 7)
}

# -----------------------------------------------------------------------------
# 4. SUMMARIZED PARAMETER MAE PLOT (parameters.png)
# -----------------------------------------------------------------------------
mae_data <- est_long %>%
  group_by(reference_method, param_set, comparison_method, sim_num) %>%
  summarise(
    state_mae = mean(abs(error[parameter %in% state_params]), na.rm = TRUE),
    obs_mae   = mean(abs(error[parameter %in% obs_params]), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(facet_labeller, by = "param_set")

# Top row: State MAE with titles and independent Y
p_mae_state <- ggplot(mae_data, aes(x = comparison_method, y = state_mae, fill = comparison_method)) +
  geom_boxplot() + 
  facet_wrap(~col_title, nrow = 1, labeller = label_parsed, scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(y = "State MAE", x = NULL)

# Bottom row: Obs MAE without titles and independent Y
p_mae_obs <- ggplot(mae_data, aes(x = comparison_method, y = obs_mae, fill = comparison_method)) +
  geom_boxplot() + 
  facet_wrap(~param_set, nrow = 1, scales = "free_y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "Obs MAE", x = "Clustering Method")

ggsave(file.path(plot_output_dir, "parameters.png"), p_mae_state / p_mae_obs, width = 12, height = 8)

# -----------------------------------------------------------------------------
# 5. SUMMARIZED PREDICTIONS PLOT (predictions.png)
# -----------------------------------------------------------------------------
pred_perf_df$comparison_method <- factor(pred_perf_df$comparison_method, levels = method_order)
pred_data <- pred_perf_df %>%
  left_join(facet_labeller, by = "param_set")

# Top row: AUC with titles and independent Y
p_auc <- ggplot(pred_data, aes(x = comparison_method, y = auc, fill = comparison_method)) +
  geom_boxplot() + 
  facet_wrap(~col_title, nrow = 1, labeller = label_parsed, scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_blank(), legend.position = "none") +
  labs(y = "AUC", x = NULL)

# Bottom row: AUPRC without titles and independent Y
p_auprc <- ggplot(pred_data, aes(x = comparison_method, y = auprc, fill = comparison_method)) +
  geom_boxplot() + 
  facet_wrap(~param_set, nrow = 1, scales = "free_y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "AUPRC", x = "Clustering Method")

ggsave(file.path(plot_output_dir, "predictions.png"), p_auc / p_auprc, width = 12, height = 8)

cat("All plots generated in:", plot_output_dir, "\n")