library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stringr)

# -------------------------------------------------------------------------
# Setup Directories
# -------------------------------------------------------------------------
input_dir <- file.path("species_experiments", "output", "points")
output_dir <- file.path("species_experiments", "output", "points", "plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------------------------------------------------------
# 1. Load Data & Preprocess
# -------------------------------------------------------------------------
data_path <- file.path(input_dir, "points_results.csv")
if (!file.exists(data_path)) {
  # Fallback if running from the output folder directly
  data_path <- "points_results.csv"
}

data <- read.csv(data_path, header = TRUE)

# --- Load species descriptions ---
species_map_path <- "./species_descr_ext.csv"
if (file.exists(species_map_path)) {
  species_map <- read.csv(species_map_path)
  data <- data %>%
    left_join(species_map %>% select(Species, Abbreviation), by = c("species" = "Abbreviation")) %>%
    mutate(species = ifelse(!is.na(Species), Species, species)) %>%
    select(-Species)
} else {
  warning("species_descr_ext.csv not found. Using abbreviations.")
}

# -------------------------------------------------------------------------
# 2. Define Color Palettes
# -------------------------------------------------------------------------

# Palette for Experiment A (3 methods)
colors_exp_a <- c(
  "lat-long" = "grey30",
  "1to10" = "#E7B800",
  "2to10" = "#FC4E07"
)

# Palette for Experiment B (18 methods)
# Distinct colors to differentiate the 18 combinations
colors_exp_b <- c(
  # 100m Group
  "occu-lat-long-100m"  = "#A6CEE3", "occuN-lat-long-100m" = "#1F78B4",
  "occu-1to10-100m"     = "#B2DF8A", "occuN-1to10-100m"    = "#33A02C",
  "occu-2to10-100m"     = "#FB9A99", "occuN-2to10-100m"    = "#E31A1C",
  
  # 200m Group
  "occu-lat-long-200m"  = "#FDBF6F", "occuN-lat-long-200m" = "#FF7F00",
  "occu-1to10-200m"     = "#CAB2D6", "occuN-1to10-200m"    = "#6A3D9A",
  "occu-2to10-200m"     = "#FFFF99", "occuN-2to10-200m"    = "#B15928",
  
  # 500m Group
  "occu-lat-long-500m"  = "#8DD3C7", "occuN-lat-long-500m" = "#BEBADA",
  "occu-1to10-500m"     = "#FB8072", "occuN-1to10-500m"    = "#80B1D3",
  "occu-2to10-500m"     = "#FDB462", "occuN-2to10-500m"    = "#B3DE69"
)

# -------------------------------------------------------------------------
# 3. Helper Functions
# -------------------------------------------------------------------------

# Function to calculate % improvement against a FIXED baseline method
calculate_improvement_fixed <- function(df, metric_col, baseline_method_name, group_cols) {
  
  # 1. Extract Baseline Rows (The Reference Standard)
  baseline_df <- df %>%
    filter(method_label == baseline_method_name) %>%
    select(species, test_repeat, baseline_val = all_of(metric_col))
  
  if(nrow(baseline_df) == 0) stop(paste("Baseline method not found:", baseline_method_name))
  
  # 2. Join Baseline to EVERYTHING (including itself)
  #    We join by species and test_repeat to ensure matched comparison
  df_diff <- df %>%
    inner_join(baseline_df, by = c("species", "test_repeat")) %>%
    mutate(perc_diff = ((.data[[metric_col]] - baseline_val) / baseline_val) * 100) %>%
    select(species, test_repeat, method_label, perc_diff)
  
  return(df_diff)
}

# Aggregation Helpers
aggregate_by_species <- function(df_diff) {
  df_diff %>%
    group_by(species, method_label) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

aggregate_by_repeat <- function(df_diff) {
  df_diff %>%
    group_by(test_repeat, method_label) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

# Plotting Function
plot_improvement <- function(df, x_col, y_col, fill_col, y_lab, palette, output_path, sort_x = TRUE) {
  
  if(sort_x) {
    # Order by mean value
    order <- df %>%
      group_by(.data[[x_col]]) %>%
      summarise(mean_val = mean(.data[[y_col]], na.rm = TRUE)) %>%
      arrange(mean_val) %>%
      pull(.data[[x_col]])
    df[[x_col]] <- factor(df[[x_col]], levels = order)
  }
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'darkgrey') +
    geom_boxplot(outlier.size = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, col = "white", bg = "darkred") +
    scale_fill_manual(values = palette) +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none" # Hide legend if x-axis labels are sufficient
    ) +
    labs(y = y_lab, x = "Method")
  
  ggsave(output_path, plot = p, width = 10, height = 7, dpi = 300)
}

# Raw Performance Plotter
plot_raw <- function(df, x_col, y_col, fill_col, y_lab, palette, output_path) {
  
  # Order species by performance
  sp_order <- df %>%
    group_by(species) %>%
    summarise(mean_val = mean(.data[[y_col]], na.rm=T)) %>%
    arrange(mean_val) %>%
    pull(species)
  df$species <- factor(df$species, levels = sp_order)
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
    geom_boxplot(outlier.size = 0.5, lwd = 0.3) +
    theme_classic() +
    coord_flip() +
    scale_fill_manual(values = palette) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    labs(x = "Species", y = y_lab)
  
  ggsave(output_path, plot = p, width = 8, height = 12, dpi = 300)
}


# =========================================================================
# PART 1: EXPERIMENT A (Unbuffered, 0m, occu only)
# =========================================================================
cat("--- Processing Experiment A (0m) ---\n")

# Filter Data
df_a <- data %>% 
  filter(buffer == 0) %>%
  mutate(method_label = method) # Use simple names: lat-long, 1to10, 2to10

# 1. Raw Plots
plot_raw(df_a, "species", "auc", "method_label", "AUC (0m)", colors_exp_a, 
         file.path(output_dir, "auc_0m.png"))
plot_raw(df_a, "species", "auprc", "method_label", "AUPRC (0m)", colors_exp_a, 
         file.path(output_dir, "auprc_0m.png"))

# 2. Percentage Difference (Baseline: lat-long)
auc_diff_a <- calculate_improvement_fixed(df_a, "auc", "lat-long")
auprc_diff_a <- calculate_improvement_fixed(df_a, "auprc", "lat-long")

# 3. Aggregations & Plots
# Species Aggregation
plot_improvement(aggregate_by_species(auc_diff_a), "method_label", "mean_perc_diff", "method_label", 
                 "% AUC Improvement (vs lat-long 0m)", colors_exp_a, file.path(output_dir, "auc_perc_diff_species_0m.png"))
plot_improvement(aggregate_by_species(auprc_diff_a), "method_label", "mean_perc_diff", "method_label", 
                 "% AUPRC Improvement (vs lat-long 0m)", colors_exp_a, file.path(output_dir, "auprc_perc_diff_species_0m.png"))

# Repeats Aggregation
plot_improvement(aggregate_by_repeat(auc_diff_a), "method_label", "mean_perc_diff", "method_label", 
                 "% AUC Improvement (vs lat-long 0m)", colors_exp_a, file.path(output_dir, "auc_perc_diff_repeats_0m.png"))
plot_improvement(aggregate_by_repeat(auprc_diff_a), "method_label", "mean_perc_diff", "method_label", 
                 "% AUPRC Improvement (vs lat-long 0m)", colors_exp_a, file.path(output_dir, "auprc_perc_diff_repeats_0m.png"))


# =========================================================================
# PART 2: EXPERIMENT B (Buffered, 18 Methods)
# =========================================================================
cat("--- Processing Experiment B (Buffered) ---\n")

# Filter & Create Label
df_b <- data %>%
  filter(buffer > 0) %>%
  mutate(
    # Construct label: "model-method-buffer"
    method_label = paste(model, method, paste0(buffer, "m"), sep = "-")
  )

# Define Baseline
baseline_b <- "occuN-lat-long-200m"

# 1. Raw Plots
plot_raw(df_b, "species", "auc", "method_label", "AUC (Buffered)", colors_exp_b, 
         file.path(output_dir, "auc_buffered.png"))
plot_raw(df_b, "species", "auprc", "method_label", "AUPRC (Buffered)", colors_exp_b, 
         file.path(output_dir, "auprc_buffered.png"))

# 2. Percentage Difference (Baseline: occuN-lat-long-200m)
auc_diff_b <- calculate_improvement_fixed(df_b, "auc", baseline_b)
auprc_diff_b <- calculate_improvement_fixed(df_b, "auprc", baseline_b)

# 3. Aggregations & Plots
# Note: X-axis sorting is enabled (sort_x=TRUE) to arrange bars by performance, 
# but you can set sort_x=FALSE to keep alphabetical/group order if preferred.

# Species Aggregation
plot_improvement(aggregate_by_species(auc_diff_b), "method_label", "mean_perc_diff", "method_label", 
                 paste("% AUC Improvement (vs", baseline_b, ")"), colors_exp_b, 
                 file.path(output_dir, "auc_perc_diff_species_buffered.png"))

plot_improvement(aggregate_by_species(auprc_diff_b), "method_label", "mean_perc_diff", "method_label", 
                 paste("% AUPRC Improvement (vs", baseline_b, ")"), colors_exp_b, 
                 file.path(output_dir, "auprc_perc_diff_species_buffered.png"))

# Repeats Aggregation
plot_improvement(aggregate_by_repeat(auc_diff_b), "method_label", "mean_perc_diff", "method_label", 
                 paste("% AUC Improvement (vs", baseline_b, ")"), colors_exp_b, 
                 file.path(output_dir, "auc_perc_diff_repeat_buffered.png")) # Matches your requested filename (repeat vs repeats)

plot_improvement(aggregate_by_repeat(auprc_diff_b), "method_label", "mean_perc_diff", "method_label", 
                 paste("% AUPRC Improvement (vs", baseline_b, ")"), colors_exp_b, 
                 file.path(output_dir, "auprc_perc_diff_repeats_buffered.png"))

cat("Done. All plots generated in:", output_dir, "\n")