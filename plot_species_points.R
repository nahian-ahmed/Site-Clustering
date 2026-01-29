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

calculate_improvement_local <- function(df, metric_col) {
  baseline_df <- df %>%
    filter(method == "lat-long") %>%
    select(species, buffer, model, test_repeat, baseline_val = all_of(metric_col))
  
  df_diff <- df %>%
    inner_join(baseline_df, by = c("species", "buffer", "model", "test_repeat")) %>%
    mutate(perc_diff = ((.data[[metric_col]] - baseline_val) / baseline_val) * 100) %>%
    select(species, buffer, model, test_repeat, method, perc_diff)
  
  return(df_diff)
}

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

df_a <- data %>% filter(buffer == 0)
colors_a <- c("lat-long" = "navy", "1to10" = "cyan", "2to10" = "pink")

# --- Raw Plots A ---
plot_raw_a <- function(df, y_col, y_lab, output_filename) {
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

# --- Perc Diff Plots A ---
plot_perc_a <- function(df, y_lab, output_filename) {
  # FORCE X-AXIS ORDER
  df$method <- factor(df$method, levels = c("lat-long", "1to10", "2to10"))
  
  p <- ggplot(df, aes(x = method, y = mean_perc_diff, fill = method)) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'darkgrey') +
    geom_boxplot() +
    # MAROON DIAMOND, WHITE BORDER
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "maroon", color = "white") +
    scale_fill_manual(values = colors_a) +
    labs(y = y_lab, x = "Method") +
    theme(legend.position = "none")
  
  # REDUCED HEIGHT (7 -> 5)
  ggsave(output_filename, plot = p, width = 6, height = 6, dpi = 300)
}

# Experiment A Plots
# Y-axis label changed to remove "(0m)"
plot_raw_a(df_a, "auc", "AUC", file.path(output_dir, "occu_0m_auc.png"))
plot_raw_a(df_a, "auprc", "AUPRC", file.path(output_dir, "occu_0m_auprc.png"))

auc_diff_a <- calculate_improvement_local(df_a, "auc")
auprc_diff_a <- calculate_improvement_local(df_a, "auprc")

# Y-axis labels updated
plot_perc_a(aggregate_by_species(auc_diff_a), "% AUC Improvement over lat-long", file.path(output_dir, "occu_0m_auc_perc_diff_species.png"))
plot_perc_a(aggregate_by_repeat(auc_diff_a), "% AUC Improvement over lat-long", file.path(output_dir, "occu_0m_auc_perc_diff_repeats.png"))
plot_perc_a(aggregate_by_species(auprc_diff_a), "% AUPRC Improvement over lat-long", file.path(output_dir, "occu_0m_auprc_perc_diff_species.png"))
plot_perc_a(aggregate_by_repeat(auprc_diff_a), "% AUPRC Improvement over lat-long", file.path(output_dir, "occu_0m_auprc_perc_diff_repeats.png"))


# -------------------------------------------------------------------------
# 4. Experiment B (Buffered)
# -------------------------------------------------------------------------
cat("--- Processing Experiment B (Buffered) ---\n")

colors_b <- c("lat-long" = "navy", "1to10" = "cyan", "2to10" = "pink")

# --- Raw Plots B ---
df_b_raw <- data %>% 
  filter(buffer > 0) %>%
  mutate(
    model = factor(model, levels = c("occuN", "occu")),
    buffer = factor(buffer, levels = c(100, 200, 500)),
    # Combined label for 1x6 grid
    panel_label = factor(paste(model, paste0(buffer, "m"), sep = "-"), 
                         levels = c("occuN-100m", "occuN-200m", "occuN-500m", 
                                    "occu-100m", "occu-200m", "occu-500m"))
  )

plot_raw_b <- function(df, y_col, y_lab, output_filename) {
  
  sp_order <- df %>% group_by(species) %>% summarise(m=mean(.data[[y_col]])) %>% arrange(m) %>% pull(species)
  df$species <- factor(df$species, levels = sp_order)
  
  p <- ggplot(df, aes(x = species, y = .data[[y_col]], fill = method)) +
    geom_boxplot(outlier.size = 0.5, lwd = 0.2) +
    theme_classic() +
    coord_flip() +
    facet_wrap(~panel_label, nrow = 1, scales = "free_x") + 
    scale_fill_manual(values = colors_b) +
    labs(x = "Species", y = y_lab) +
    theme(
      legend.position = "bottom", 
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  ggsave(output_filename, plot = p, width = 18, height = 12, dpi = 300)
}

plot_raw_b(df_b_raw, "auc", "AUC (Buffered)", file.path(output_dir, "buffered_auc.png"))
plot_raw_b(df_b_raw, "auprc", "AUPRC (Buffered)", file.path(output_dir, "buffered_auprc.png"))


# --- Perc Diff Plots B ---
df_b_diff <- data %>% filter(buffer > 0)

auc_diff_b <- calculate_improvement_local(df_b_diff, "auc")
auprc_diff_b <- calculate_improvement_local(df_b_diff, "auprc")

# Helper to format the buffer column for plotting
format_buffer_label <- function(df) {
  df %>% mutate(
    model = factor(model, levels = c("occuN", "occu")),
    buffer_label = factor(paste0(buffer, "m"), levels = c("100m", "200m", "500m"))
  )
}

plot_perc_b <- function(df, y_lab, output_filename) {
  
  df_formatted <- format_buffer_label(df)
  
  # FORCE X-AXIS ORDER
  df_formatted$method <- factor(df_formatted$method, levels = c("lat-long", "1to10", "2to10"))
  
  p <- ggplot(df_formatted, aes(x = method, y = mean_perc_diff, fill = method)) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'darkgrey') +
    geom_boxplot() +
    # MAROON DIAMOND, WHITE BORDER
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "maroon", color = "white") +
    facet_grid(model ~ buffer_label, scales = "fixed") + 
    scale_fill_manual(values = colors_b) +
    labs(y = y_lab, x = "Method") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  # REDUCED HEIGHT (8 -> 6)
  ggsave(output_filename, plot = p, width = 10, height = 7, dpi = 300)
}

agg_auc_species <- aggregate_by_species(auc_diff_b)
agg_auc_repeat <- aggregate_by_repeat(auc_diff_b)
agg_auprc_species <- aggregate_by_species(auprc_diff_b)
agg_auprc_repeat <- aggregate_by_repeat(auprc_diff_b)

# Y-axis labels updated
plot_perc_b(agg_auc_species, "% AUC Improvement over lat-long", file.path(output_dir, "buffered_auc_perc_diff_species.png"))
plot_perc_b(agg_auc_repeat, "% AUC Improvement over lat-long", file.path(output_dir, "buffered_auc_perc_diff_repeats.png"))
plot_perc_b(agg_auprc_species, "% AUPRC Improvement over lat-long", file.path(output_dir, "buffered_auprc_perc_diff_species.png"))
plot_perc_b(agg_auprc_repeat, "% AUPRC Improvement over lat-long", file.path(output_dir, "buffered_auprc_perc_diff_repeats.png"))

cat("Done. All plots generated in:", output_dir, "\n")