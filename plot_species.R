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

# --- Load species descriptions and replace abbreviations with full names ---
species_map <- read.csv("./species_descr_ext.csv")

# Join to get full names. 
# We match data$species (abbreviations) with species_map$Abbreviation
data <- data %>%
  left_join(species_map %>% select(Species, Abbreviation), by = c("species" = "Abbreviation")) %>%
  # If a match is found, replace 'species' with the full name 'Species'
  mutate(species = ifelse(!is.na(Species), Species, species)) %>%
  select(-Species)
# -------------------------------------------------------------------------

# Define the manual color scale
colors <- c("forestgreen", "darkgrey", "red", "blue", "yellow", "orange", 
            "green", "pink", "cyan", "navy", "brown")

names(colors) <- c("BayesOptClustGeo", "DBSC", "best-clustGeo", "1-per-UL", "SVS", "rounded-4",
                   "1-kmSq", "2to10-sameObs", "2to10", "1to10", "lat-long")

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
final_data_auc <- bind_rows(other_data, best_cg_data_auc)
final_data_auprc <- bind_rows(other_data, best_cg_data_auprc)

# -------------------------------------------------------------------------
# (2) Plot Raw AUC and AUPRC
# -------------------------------------------------------------------------

plot_raw_performance <- function(df, metric_col, y_label, output_filename) {
  
  species_order <- df %>%
    group_by(species) %>%
    summarise(metric_mean = mean(.data[[metric_col]], na.rm = TRUE)) %>%
    arrange(metric_mean) %>%
    pull(species)
  
  df$species <- factor(df$species, levels = species_order)
  
  alg_order <- c("2to10", "2to10-sameObs", "1to10", "1-kmSq", "lat-long", 
                 "rounded-4", "SVS", "1-per-UL", "best-clustGeo", "DBSC", "BayesOptClustGeo")
  df <- df %>% filter(method %in% alg_order)
  df$method <- factor(df$method, levels = alg_order)
  
  p <- ggplot(df, aes(x = species, y = .data[[metric_col]], fill = method)) +
    geom_boxplot(outlier.size = 0.5, lwd = 0.3) +
    theme_classic() +
    coord_flip()+
    scale_fill_manual(values = colors) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    labs(x = "Species", y = y_label)
  
  ggsave(output_filename, plot = p, width = 7, height = 12, dpi = 300)
}

plot_raw_performance(final_data_auc, "auc", "AUC", 
                     file.path(output_plot_dir, "auc.png"))

plot_raw_performance(final_data_auprc, "auprc", "AUPRC", 
                     file.path(output_plot_dir, "auprc.png"))


# -------------------------------------------------------------------------
# (3) Calculate Percentage Improvement
# -------------------------------------------------------------------------

# Helper to prepare the diff dataframe
prepare_diff_data <- function(df, metric_col) {
  # Pivot to wide
  df_wide <- df %>%
    select(species, test_repeat, method, all_of(metric_col)) %>%
    tidyr::pivot_wider(names_from = method, values_from = all_of(metric_col))
  
  if(!"lat-long" %in% colnames(df_wide)) stop("lat-long method missing from data")
  
  # Store baseline 
  baseline_values <- df_wide[["lat-long"]]
  
  methods_to_compare <- setdiff(names(df_wide), c("species", "test_repeat"))
  
  df_diff <- df_wide
  for(m in methods_to_compare) {
    # Calculate % diff
    df_diff[[m]] <- ((df_diff[[m]] - baseline_values) / baseline_values) * 100
  }
  
  # Pivot back to long
  df_diff_long <- df_diff %>%
    tidyr::pivot_longer(cols = all_of(methods_to_compare), 
                        names_to = "method", 
                        values_to = "perc_diff")
  
  return(df_diff_long)
}

# --- Aggregation Strategy 1: Average 25 repeats per Species (31 points per method) ---
aggregate_by_species <- function(df_long) {
  df_long %>%
    group_by(species, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

# --- Aggregation Strategy 2: Average 31 species per Repeat (25 points per method) ---
aggregate_by_repeat <- function(df_long) {
  df_long %>%
    group_by(test_repeat, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

# Prepare base diff data
auc_diff_raw   <- prepare_diff_data(final_data_auc, "auc")
auprc_diff_raw <- prepare_diff_data(final_data_auprc, "auprc")

# Create the two aggregated datasets for AUC
auc_by_species <- aggregate_by_species(auc_diff_raw)
auc_by_repeat  <- aggregate_by_repeat(auc_diff_raw)

# Create the two aggregated datasets for AUPRC
auprc_by_species <- aggregate_by_species(auprc_diff_raw)
auprc_by_repeat  <- aggregate_by_repeat(auprc_diff_raw)

# -------------------------------------------------------------------------
# (4) Plot Percentage Improvement
# -------------------------------------------------------------------------

plot_improvement <- function(df, y_label, output_filename) {
  
  # Order methods by median improvement
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
  
  ggsave(output_filename, plot = p, width = 6, height = 7, dpi = 300)
}

# --- PLOT SET 1: Distribution over Species (averaged over repeats) ---
plot_improvement(auc_by_species, 
                 "Average % AUC Improvement (Average over 25 Repeats per Species)",
                 file.path(output_plot_dir, "auc_perc_diff_species.png"))

plot_improvement(auprc_by_species, 
                 "Average % AUPRC Improvement (Average over 25 Repeats per Species)",
                 file.path(output_plot_dir, "auprc_perc_diff_species.png"))

# --- PLOT SET 2: Distribution over Repeats (averaged over species) ---
plot_improvement(auc_by_repeat, 
                 "Average % AUC Improvement (Average over 31 Species per Repeat)",
                 file.path(output_plot_dir, "auc_perc_diff_repeats.png"))

plot_improvement(auprc_by_repeat, 
                 "Average % AUPRC Improvement (Average over 31 Species per Repeat)",
                 file.path(output_plot_dir, "auprc_perc_diff_repeats.png"))

print("Plots generated successfully in simulation_experiments/output/plots/")