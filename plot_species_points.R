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
# 1. Load Data & Preprocess
# -------------------------------------------------------------------------
data_path <- file.path(input_dir, "points_results.csv")
if (!file.exists(data_path)) {
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
# 2. Helper Functions
# -------------------------------------------------------------------------

# Calculate % Improvement (Local Baseline)
# Computes diff relative to 'lat-long' within the SAME buffer/model/repeat group
calculate_improvement_local <- function(df, metric_col) {
  
  # 1. Extract Baseline (lat-long) for every group
  baseline_df <- df %>%
    filter(method == "lat-long") %>%
    select(species, buffer, model, test_repeat, baseline_val = all_of(metric_col))
  
  # 2. Join Baseline back to the full dataset
  #    This ensures every row (including lat-long itself) gets a baseline_val
  df_diff <- df %>%
    inner_join(baseline_df, by = c("species", "buffer", "model", "test_repeat")) %>%
    mutate(perc_diff = ((.data[[metric_col]] - baseline_val) / baseline_val) * 100) %>%
    select(species, buffer, model, test_repeat, method, perc_diff)
  
  return(df_diff)
}

# Aggregation Helpers
aggregate_by_species <- function(df_diff) {
  df_diff %>%
    group_by(species, buffer, model, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

aggregate_by_repeat <- function(df_diff) {
  df_diff %>%
    group_by(test_repeat, buffer, model, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

# -------------------------------------------------------------------------
# 3. Experiment A (0m, Unbuffered)
# -------------------------------------------------------------------------
cat("--- Processing Experiment A (0m) ---\n")

# Filter Data
df_a <- data %>% filter(buffer == 0)

# Colors for A
colors_a <- c("lat-long" = "navy", "1to10" = "cyan", "2to10" = "pink")

# --- Plotting Function for A ---
plot_raw_a <- function(df, y_col, y_lab, output_filename) {
  # Sort Species
  sp_order <- df %>% group_by(species) %>% summarise(m=mean(.data[[y_col]])) %>% arrange(m) %>% pull(species)
  df$species <- factor(df$species, levels = sp_order)
  
  p <- ggplot(df, aes(x = species, y = .data[[y_col]], fill = method)) +
    geom_boxplot(outlier.size = 0.5, lwd = 0.3) +
    theme_classic() +
    coord_flip() +
    scale_fill_manual(values = colors_a) +
    labs(x = "Species", y = y_lab) +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(output_filename, plot = p, width = 8, height = 12, dpi = 300)
}

plot_perc_a <- function(df, y_lab, output_filename) {
  # Sort Method
  m_order <- df %>% group_by(method) %>% summarise(m=mean(mean_perc_diff)) %>% arrange(m) %>% pull(method)
  df$method <- factor(df$method, levels = m_order)
  
  p <- ggplot(df, aes(x = method, y = mean_perc_diff, fill = method)) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    scale_fill_manual(values = colors_a) +
    labs(y = y_lab, x = "Method") +
    theme(legend.position = "none")
  
  ggsave(output_filename, plot = p, width = 6, height = 7, dpi = 300)
}

# 1. Raw Plots
plot_raw_a(df_a, "auc", "AUC (0m)", file.path(output_dir, "occu_0m_auc.png"))
plot_raw_a(df_a, "auprc", "AUPRC (0m)", file.path(output_dir, "occu_0m_auprc.png"))

# 2. Perc Diff
auc_diff_a <- calculate_improvement_local(df_a, "auc")
auprc_diff_a <- calculate_improvement_local(df_a, "auprc")

plot_perc_a(aggregate_by_species(auc_diff_a), "% AUC Improvement (vs lat-long)", file.path(output_dir, "occu_0m_auc_perc_diff_species.png"))
plot_perc_a(aggregate_by_repeat(auc_diff_a), "% AUC Improvement (vs lat-long)", file.path(output_dir, "occu_0m_auc_perc_diff_repeats.png"))
plot_perc_a(aggregate_by_species(auprc_diff_a), "% AUPRC Improvement (vs lat-long)", file.path(output_dir, "occu_0m_auprc_perc_diff_species.png"))
plot_perc_a(aggregate_by_repeat(auprc_diff_a), "% AUPRC Improvement (vs lat-long)", file.path(output_dir, "occu_0m_auprc_perc_diff_repeats.png"))


# -------------------------------------------------------------------------
# 4. Experiment B (Buffered)
# -------------------------------------------------------------------------
cat("--- Processing Experiment B (Buffered) ---\n")

# Filter Data & Factor Ordering
df_b <- data %>% 
  filter(buffer > 0) %>%
  mutate(
    # Ensure Model factor level order for Rows (occuN top, occu bottom)
    model = factor(model, levels = c("occuN", "occu")),
    # Ensure Buffer factor level order for Cols (100, 200, 500)
    buffer = factor(buffer, levels = c(100, 200, 500)),
    # Create combined label for 1x6 plot ordering
    panel_label = factor(paste(model, buffer), 
                         levels = c("occuN 100", "occuN 200", "occuN 500", 
                                    "occu 100", "occu 200", "occu 500"))
  )

# Colors for B
colors_b <- c("lat-long" = "navy", "1to10" = "cyan", "2to10" = "pink")

# --- Plotting Function for B (Raw: 1x6 Grid) ---
plot_raw_b <- function(df, y_col, y_lab, output_filename) {
  
  # Sort Species by overall mean
  sp_order <- df %>% group_by(species) %>% summarise(m=mean(.data[[y_col]])) %>% arrange(m) %>% pull(species)
  df$species <- factor(df$species, levels = sp_order)
  
  p <- ggplot(df, aes(x = species, y = .data[[y_col]], fill = method)) +
    geom_boxplot(outlier.size = 0.5, lwd = 0.2) +
    theme_classic() +
    coord_flip() +
    # 1 Row, 6 Columns Layout
    facet_wrap(~panel_label, nrow = 1, scales = "free_x") + 
    scale_fill_manual(values = colors_b) +
    labs(x = "Species", y = y_lab) +
    theme(
      legend.position = "bottom", 
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  # Width increased to accommodate 6 columns
  ggsave(output_filename, plot = p, width = 18, height = 12, dpi = 300)
}

# --- Plotting Function for B (Perc Diff: 2x3 Grid) ---
plot_perc_b <- function(df, y_lab, output_filename) {
  
  # We want lat-long included (it will be 0)
  
  p <- ggplot(df, aes(x = method, y = mean_perc_diff, fill = method)) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
    # 2 Rows (Model) x 3 Columns (Buffer)
    facet_grid(model ~ buffer, scales = "fixed") + 
    scale_fill_manual(values = colors_b) +
    labs(y = y_lab, x = "Method") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(output_filename, plot = p, width = 10, height = 8, dpi = 300)
}

# 1. Raw Plots
plot_raw_b(df_b, "auc", "AUC (Buffered)", file.path(output_dir, "buffered_auc.png"))
plot_raw_b(df_b, "auprc", "AUPRC (Buffered)", file.path(output_dir, "buffered_auprc.png"))

# 2. Perc Diff
auc_diff_b <- calculate_improvement_local(df_b, "auc")
auprc_diff_b <- calculate_improvement_local(df_b, "auprc")

# Ensure aggregation preserves the grouping factors for faceting
agg_auc_species <- aggregate_by_species(auc_diff_b)
agg_auc_repeat <- aggregate_by_repeat(auc_diff_b)
agg_auprc_species <- aggregate_by_species(auprc_diff_b)
agg_auprc_repeat <- aggregate_by_repeat(auprc_diff_b)

plot_perc_b(agg_auc_species, "% AUC Improvement", file.path(output_dir, "buffered_auc_perc_diff_species.png"))
plot_perc_b(agg_auc_repeat, "% AUC Improvement", file.path(output_dir, "buffered_auc_perc_diff_repeats.png"))
plot_perc_b(agg_auprc_species, "% AUPRC Improvement", file.path(output_dir, "buffered_auprc_perc_diff_species.png"))
plot_perc_b(agg_auprc_repeat, "% AUPRC Improvement", file.path(output_dir, "buffered_auprc_perc_diff_repeats.png"))

cat("Done. All plots generated in:", output_dir, "\n")