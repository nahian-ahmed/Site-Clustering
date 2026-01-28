library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# -------------------------------------------------------------------------
# Setup Directories
# -------------------------------------------------------------------------
input_dir <- file.path("species_experiments", "output", "points")
output_dir <- file.path("species_experiments", "output", "points", "plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------------------------------------------------------
# 1. Load Data
# -------------------------------------------------------------------------
data_path <- file.path(input_dir, "points_results.csv")
if (!file.exists(data_path)) {
  # Fallback if running from the output folder directly
  data_path <- "points_results.csv"
}

data <- read.csv(data_path, header = TRUE)

# --- Load species descriptions and replace abbreviations with full names ---
species_map_path <- "./species_descr_ext.csv"
if (file.exists(species_map_path)) {
  species_map <- read.csv(species_map_path)
  
  # Join to get full names
  data <- data %>%
    left_join(species_map %>% select(Species, Abbreviation), by = c("species" = "Abbreviation")) %>%
    mutate(species = ifelse(!is.na(Species), Species, species)) %>%
    select(-Species)
} else {
  warning("species_descr_ext.csv not found. Using abbreviations.")
}

# -------------------------------------------------------------------------
# 2. Define Colors
# -------------------------------------------------------------------------
colors <- c("forestgreen", "darkgrey", "red", "blue", "yellow", "orange", 
            "green", "pink", "cyan", "navy", "brown")

names(colors) <- c("BayesOptClustGeo", "DBSC", "best-clustGeo", "best-SLIC", "rounded-4",
                   "1-kmSq", "2to10-sameObs", "2to10", "1to10", "lat-long")

# -------------------------------------------------------------------------
# 3. Calculate Percentage Improvement (Robust Join Method)
# -------------------------------------------------------------------------

calculate_improvement <- function(df, metric_col) {
  
  # 1. Extract the Baseline (lat-long) for every condition
  #    We group by Species, Buffer, Model, and Repeat to ensure apples-to-apples comparison.
  baseline_df <- df %>%
    filter(method == "lat-long") %>%
    select(species, buffer, model, test_repeat, baseline_val = all_of(metric_col))
  
  # 2. Extract the Methods to Compare (excluding lat-long)
  compare_df <- df %>%
    filter(method != "lat-long") %>%
    select(species, buffer, model, test_repeat, method, val = all_of(metric_col))
  
  # 3. Join and Calculate Diff
  #    This automatically matches '2to10' (occuN) with 'lat-long' (occuN)
  #    and '2to10' (occu) with 'lat-long' (occu)
  result_df <- compare_df %>%
    inner_join(baseline_df, by = c("species", "buffer", "model", "test_repeat")) %>%
    mutate(perc_diff = ((val - baseline_val) / baseline_val) * 100) %>%
    select(species, buffer, model, test_repeat, method, perc_diff)
  
  return(result_df)
}

# --- Aggregation Strategy 1: Average over Repeats (Group by Species) ---
aggregate_by_species <- function(df_diff) {
  df_diff %>%
    # Average across repeats, buffers, and models to get one summary value per Species per Method
    group_by(species, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

# --- Aggregation Strategy 2: Average over Species (Group by Repeat) ---
aggregate_by_repeat <- function(df_diff) {
  df_diff %>%
    # Average across species, buffers, and models to get one summary value per Repeat per Method
    group_by(test_repeat, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

# 1. Compute raw differences
auc_diff_raw   <- calculate_improvement(data, "auc")
auprc_diff_raw <- calculate_improvement(data, "auprc")

# 2. Create the aggregated datasets for AUC
auc_by_species <- aggregate_by_species(auc_diff_raw)
auc_by_repeat  <- aggregate_by_repeat(auc_diff_raw)

# 3. Create the aggregated datasets for AUPRC
auprc_by_species <- aggregate_by_species(auprc_diff_raw)
auprc_by_repeat  <- aggregate_by_repeat(auprc_diff_raw)

# -------------------------------------------------------------------------
# 4. Plot Percentage Improvement
# -------------------------------------------------------------------------

plot_improvement <- function(df, y_label, output_filename) {
  
  # Determine order by mean improvement
  method_order <- df %>%
    group_by(method) %>%
    summarise(mean_val = mean(mean_perc_diff, na.rm = TRUE)) %>%
    arrange(mean_val) %>%
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
  
  ggsave(output_filename, plot = p, width = 6, height = 7, dpi = 300)
}

# --- PLOT SET 1: Distribution over Species ---
plot_improvement(auc_by_species, 
                 "Average % AUC Improvement (Avged over Repeats & Models)",
                 file.path(output_dir, "auc_perc_diff_species.png"))

plot_improvement(auprc_by_species, 
                 "Average % AUPRC Improvement (Avged over Repeats & Models)",
                 file.path(output_dir, "auprc_perc_diff_species.png"))

# --- PLOT SET 2: Distribution over Repeats ---
plot_improvement(auc_by_repeat, 
                 "Average % AUC Improvement (Avged over Species & Models)",
                 file.path(output_dir, "auc_perc_diff_repeats.png"))

plot_improvement(auprc_by_repeat, 
                 "Average % AUPRC Improvement (Avged over Species & Models)",
                 file.path(output_dir, "auprc_perc_diff_repeats.png"))

print(paste("Plots generated successfully in", output_dir))