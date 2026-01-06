library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# Create output directory
# dir.create("simulation_experiments/output/plots", recursive = TRUE, showWarnings = FALSE)

output_dir <- file.path("species_experiments", "output")
output_plot_dir <- file.path(output_dir, "plots")
if (!dir.exists(output_plot_dir)) dir.create(output_plot_dir, recursive = TRUE)


# 1. Load Data
data <- read.delim(file.path(output_dir, "predictive_performance.csv"), sep = ",", header = T)

# Ensure test_repeat is a factor if needed, though for aggregation strictly it's numeric
# The unique methods are:
# "1to10", "2to10", "2to10-sameObs", "1-kmSq", "lat-long", "rounded-4", 
# "SVS", "1-per-UL", "DBSC", "BayesOptClustGeo", and various "clustGeo-X-Y"

# Define the manual color scale from your old script
# Note: I added "best-clustGeo" to the colors mapping
colors <- c("forestgreen", "darkgrey", "red", "blue", "yellow", "orange", 
            "purple", "green", "brown", "pink", "cyan", "black")
names(colors) <- c("2to10", "2to10-sameObs", "1to10", "1-kmSq", "lat-long", 
                   "rounded-4", "SVS", "1-per-UL", "DBSC", "BayesOptClustGeo", 
                   "clustGeo", "best-clustGeo")

# -------------------------------------------------------------------------
# (1) Define best-clustGeo
# -------------------------------------------------------------------------

# Separate clustGeo variants (excluding BayesOpt) from the rest
cg_data <- data %>% 
  filter(grepl("^clustGeo", method) & method != "BayesOptClustGeo")

other_data <- data %>% 
  filter(!grepl("^clustGeo", method) | method == "BayesOptClustGeo")

# Calculate mean performance over repeats for each clustGeo variant per species
cg_summary <- cg_data %>%
  group_by(species, method) %>%
  summarise(
    mean_auc = mean(auc, na.rm = TRUE),
    mean_auprc = mean(auprc, na.rm = TRUE),
    .groups = "drop"
  )

# Pick best clustGeo for AUC
best_cg_auc_methods <- cg_summary %>%
  group_by(species) %>%
  slice_max(mean_auc, n = 1, with_ties = FALSE) %>%
  select(species, best_method_auc = method)

# Pick best clustGeo for AUPRC
best_cg_auprc_methods <- cg_summary %>%
  group_by(species) %>%
  slice_max(mean_auprc, n = 1, with_ties = FALSE) %>%
  select(species, best_method_auprc = method)

# Create "best-clustGeo" rows for AUC dataset
best_cg_data_auc <- cg_data %>%
  inner_join(best_cg_auc_methods, by = c("species")) %>%
  filter(method == best_method_auc) %>%
  mutate(method = "best-clustGeo") %>%
  select(-best_method_auc)

# Create "best-clustGeo" rows for AUPRC dataset
best_cg_data_auprc <- cg_data %>%
  inner_join(best_cg_auprc_methods, by = c("species")) %>%
  filter(method == best_method_auprc) %>%
  mutate(method = "best-clustGeo") %>%
  select(-best_method_auprc)

# Combine with other data
# We need two master datasets because best-clustGeo differs for AUC vs AUPRC
final_data_auc <- bind_rows(other_data, best_cg_data_auc)
final_data_auprc <- bind_rows(other_data, best_cg_data_auprc)

# -------------------------------------------------------------------------
# (2) Plot Raw AUC and AUPRC over 25 test repeats
# -------------------------------------------------------------------------

# Helper function to plot raw performance
plot_raw_performance <- function(df, metric_col, y_label, output_filename) {
  
  # Calculate species ordering by mean median performance
  species_order <- df %>%
    group_by(species) %>%
    summarise(metric_mean = mean(.data[[metric_col]], na.rm = TRUE)) %>%
    arrange(desc(metric_mean)) %>%
    pull(species)
  
  df$species <- factor(df$species, levels = species_order)
  
  # Standard algorithm ordering
  alg_order <- c("2to10", "2to10-sameObs", "1to10", "1-kmSq", "lat-long", 
                 "rounded-4", "SVS", "1-per-UL", "best-clustGeo", "DBSC", "BayesOptClustGeo")
  # Filter to only keep methods in our list to avoid color errors if new ones exist
  df <- df %>% filter(method %in% alg_order)
  df$method <- factor(df$method, levels = alg_order)
  
  p <- ggplot(df, aes(x = species, y = .data[[metric_col]], fill = method)) +
    geom_boxplot(outlier.size = 0.5, lwd = 0.3) +
    theme_classic() +
    scale_fill_manual(values = colors) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    labs(x = "Species", y = y_label)
  
  ggsave(output_filename, plot = p, width = 14, height = 8, dpi = 400)
}

# Plot AUC
plot_raw_performance(final_data_auc, "auc", "AUC", 
                     file.path(output_plot_dir, "auc.png"))

# Plot AUPRC
plot_raw_performance(final_data_auprc, "auprc", "AUPRC", 
                     file.path(output_plot_dir, "auprc.png"))


# -------------------------------------------------------------------------
# (3) Calculate Percentage Improvement over lat-long
# -------------------------------------------------------------------------

calculate_improvement <- function(df, metric_col) {
  # Pivot to have methods as columns to easily subtract lat-long
  df_wide <- df %>%
    select(species, test_repeat, method, all_of(metric_col)) %>%
    tidyr::pivot_wider(names_from = method, values_from = all_of(metric_col))
  
  # Check if lat-long exists
  if(!"lat-long" %in% colnames(df_wide)) stop("lat-long method missing from data")
  
  # Calculate % diff for each method vs lat-long
  # Formula: ((Method - LatLong) / LatLong) * 100
  methods_to_compare <- setdiff(names(df_wide), c("species", "test_repeat", "lat-long"))
  
  df_diff <- df_wide
  for(m in methods_to_compare) {
    df_diff[[m]] <- ((df_diff[[m]] - df_diff[["lat-long"]]) / df_diff[["lat-long"]]) * 100
  }
  
  # Remove lat-long column (diff is 0) and pivot back to long
  df_diff_long <- df_diff %>%
    select(-`lat-long`) %>%
    tidyr::pivot_longer(cols = all_of(methods_to_compare), 
                        names_to = "method", 
                        values_to = "perc_diff")
  
  # AVERAGE OVER REPEATS (The "Right Way")
  # Result: One value per species per method
  df_species_avg <- df_diff_long %>%
    group_by(species, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
  
  return(df_species_avg)
}

auc_diff_df <- calculate_improvement(final_data_auc, "auc")
auprc_diff_df <- calculate_improvement(final_data_auprc, "auprc")

# -------------------------------------------------------------------------
# (4) Plot Percentage Improvement
# -------------------------------------------------------------------------

plot_improvement <- function(df, y_label, output_filename) {
  
  # Order methods by median improvement across species
  method_order <- df %>%
    group_by(method) %>%
    summarise(median_diff = median(mean_perc_diff, na.rm = TRUE)) %>%
    arrange(median_diff) %>%
    pull(method)
  
  df$method <- factor(df$method, levels = method_order)
  
  p <- ggplot(df, aes(x = method, y = mean_perc_diff, fill = method)) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, col = "white", bg = "darkred") +
    scale_fill_manual(values = colors) +
    theme(
      axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size = 11),
      axis.title = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(
      y = y_label,
      x = "Algorithm"
    )
  
  ggsave(output_filename, plot = p, width = 10, height = 8, dpi = 400)
}

plot_improvement(auc_diff_df, 
                 paste0("Average % AUC Improvement (over ", length(unique(auc_diff_df$species)), " Species)"),
                 file.path(output_plot_dir, "auc_perc_diff.png"))

plot_improvement(auprc_diff_df, 
                 paste0("Average % AUPRC Improvement (over ", length(unique(auprc_diff_df$species)), " Species)"),
                 file.path(output_plot_dir, "auprc_perc_diff.png"))

print("Plots generated successfully in simulation_experiments/output/plots/")