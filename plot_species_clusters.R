################################################################
# Plot Species Cluster-Based Experiments
# Updated: February 21, 2026
#
# INCLUDES:
# 1. Raw Performance Plots
# 2. Percentage Improvement Plots
# 3. Significance Heatmaps (% Impr, Lower Triangle)
# 4. Trait-based Analysis (Mixed-Effects on Raw AUC)
# 5. Kappa Selection Analysis & Saving Params
# 6. Maps
################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(patchwork)
library(FSA)      # For Dunn's Test
library(stringr)  # For string manipulation
library(terra)    # Loaded early to avoid masking issues
library(sf)
library(lme4)     # Mixed-effects models
library(sjPlot)   # Plotting estimates

# Create output directory
output_dir <- file.path("output", "species_experiments", "clusters")
output_plot_dir <- file.path(output_dir, "plots")
map_output_dir <- file.path(output_plot_dir, "maps")
if (!dir.exists(map_output_dir)) dir.create(map_output_dir, recursive = TRUE)

# -------------------------------------------------------------------------
# 1. Load Data
# -------------------------------------------------------------------------
data <- read.delim(file.path(output_dir, "predictive_performance.csv"), sep = ",", header = T)

# --- Load species descriptions and replace abbreviations with full names ---
species_map <- read.csv("./species_descr_ext.csv")

# Join to get full names. 
data <- data %>%
  left_join(species_map %>% select(Species, Abbreviation), by = c("species" = "Abbreviation")) %>%
  mutate(species = ifelse(!is.na(Species), Species, species)) %>%
  select(-Species)

# Define the manual color scale
colors <- c("forestgreen", "darkgrey", "red", "yellow", "orange", 
            "green", "pink", "cyan", "navy")

names(colors) <- c("BayesOptClustGeo", "DBSC", "best-clustGeo", "rounded-4",
                   "1-kmSq", "2to10-sameObs", "2to10", "1to10", "lat-long")

# -------------------------------------------------------------------------
# (1) Define best-clustGeo AND best-SLIC
# -------------------------------------------------------------------------

# 1. Separate clustGeo variants (excluding BayesOpt)
cg_data <- data %>% 
  filter(grepl("^clustGeo", method) & method != "BayesOptClustGeo")

# 3. Get the "other" data (everything else)
other_data <- data %>% 
  filter((!grepl("^clustGeo", method) | method == "BayesOptClustGeo"))

# --- Process ClustGeo ---
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

# --- Combine all data ---
final_data_auc <- bind_rows(other_data, best_cg_data_auc)
final_data_auprc <- bind_rows(other_data, best_cg_data_auprc)

# -------------------------------------------------------------------------
# (1.5) Save Best clustGeo Parameters to CSV
# -------------------------------------------------------------------------
best_params_auc <- best_cg_auc_methods %>%
  mutate(metric = "AUC") %>%
  rename(method = best_method_auc)

best_params_auprc <- best_cg_auprc_methods %>%
  mutate(metric = "AUPRC") %>%
  rename(method = best_method_auprc)

best_params_combined <- bind_rows(best_params_auc, best_params_auprc) %>%
  mutate(params_str = sub("^clustGeo-", "", method)) %>%
  separate(params_str, c("rho", "kappa"), sep = "-") %>%
  select(species, metric, method, rho, kappa)

write.csv(best_params_combined, file.path(output_plot_dir, "best-clustGeo_params.csv"), row.names = FALSE)
cat("Best clustGeo parameters saved.\n")

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
  
  # UPDATED ALGORITHM ORDER
  alg_order <- c("lat-long", "1-kmSq", "DBSC", "1to10", "2to10-sameObs", 
                 "2to10", "rounded-4", "BayesOptClustGeo", "best-clustGeo")
  
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

plot_raw_performance(final_data_auc, "auc", "AUC", file.path(output_plot_dir, "auc.png"))
plot_raw_performance(final_data_auprc, "auprc", "AUPRC", file.path(output_plot_dir, "auprc.png"))

# -------------------------------------------------------------------------
# (3) Calculate Percentage Improvement
# -------------------------------------------------------------------------
prepare_diff_data <- function(df, metric_col) {
  df_wide <- df %>%
    select(species, test_repeat, method, all_of(metric_col)) %>%
    tidyr::pivot_wider(names_from = method, values_from = all_of(metric_col))
  
  if(!"lat-long" %in% colnames(df_wide)) stop("lat-long method missing from data")
  
  baseline_values <- df_wide[["lat-long"]]
  methods_to_compare <- setdiff(names(df_wide), c("species", "test_repeat"))
  
  df_diff <- df_wide
  for(m in methods_to_compare) {
    df_diff[[m]] <- ((df_diff[[m]] - baseline_values) / baseline_values) * 100
  }
  
  df_diff_long <- df_diff %>%
    tidyr::pivot_longer(cols = all_of(methods_to_compare), 
                        names_to = "method", 
                        values_to = "perc_diff")
  return(df_diff_long)
}

aggregate_by_species <- function(df_long) {
  df_long %>%
    group_by(species, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

aggregate_by_repeat <- function(df_long) {
  df_long %>%
    group_by(test_repeat, method) %>%
    summarise(mean_perc_diff = mean(perc_diff, na.rm = TRUE), .groups = "drop")
}

auc_diff_raw   <- prepare_diff_data(final_data_auc, "auc")
auprc_diff_raw <- prepare_diff_data(final_data_auprc, "auprc")

auc_by_species <- aggregate_by_species(auc_diff_raw)
auc_by_repeat  <- aggregate_by_repeat(auc_diff_raw)
auprc_by_species <- aggregate_by_species(auprc_diff_raw)
auprc_by_repeat  <- aggregate_by_repeat(auprc_diff_raw)

# -------------------------------------------------------------------------
# (4) Plot Percentage Improvement
# -------------------------------------------------------------------------
plot_improvement <- function(df, y_label, output_filename) {
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
    labs(y = y_label, x = "Algorithm")
  
  ggsave(output_filename, plot = p, width = 6, height = 7, dpi = 300)
}

plot_improvement(auc_by_species, "Average % AUC Improvement", file.path(output_plot_dir, "auc_perc_diff_species.png"))
plot_improvement(auprc_by_species, "Average % AUPRC Improvement", file.path(output_plot_dir, "auprc_perc_diff_species.png"))
plot_improvement(auc_by_repeat, "Average % AUC Improvement", file.path(output_plot_dir, "auc_perc_diff_repeats.png"))
plot_improvement(auprc_by_repeat, "Average % AUPRC Improvement", file.path(output_plot_dir, "auprc_perc_diff_repeats.png"))

# -------------------------------------------------------------------------
# (6) Statistical Significance Heatmap
# -------------------------------------------------------------------------
cat("\n###############################################\n")
cat("GENERATING STATISTICAL SIGNIFICANCE HEATMAP (% IMPROVEMENT)\n")
cat("###############################################\n")

# Use Percentage Improvement data for the test
stats_df <- auc_by_species

# Ensure factors
stats_df$method <- factor(stats_df$method)

# Run Dunn's Test
dunn_res <- dunnTest(mean_perc_diff ~ method, data = stats_df, method = "bh")$res

# Parse comparisons (e.g., "A - B")
dunn_res <- dunn_res %>%
  separate(Comparison, into = c("Method1", "Method2"), sep = " - ")

# --- CREATE LOWER TRIANGLE ONLY ---
# 1. Create both directions to ensure full coverage
dunn_res_inv <- dunn_res
dunn_res_inv$Method1 <- dunn_res$Method2
dunn_res_inv$Method2 <- dunn_res$Method1
dunn_res_full <- bind_rows(dunn_res, dunn_res_inv)

# 2. Define performance order for axes
method_perf_order <- stats_df %>%
  group_by(method) %>%
  summarise(mean_val = mean(mean_perc_diff)) %>%
  arrange(mean_val) %>%
  pull(method)

dunn_res_full$Method1 <- factor(dunn_res_full$Method1, levels = method_perf_order)
dunn_res_full$Method2 <- factor(dunn_res_full$Method2, levels = method_perf_order)

# 3. Filter for Lower Triangle: Method2 (Y) < Method1 (X)
dunn_res_final <- dunn_res_full %>%
  filter(as.numeric(Method2) < as.numeric(Method1))
# ----------------------------------

# Add significance stars
dunn_res_final <- dunn_res_final %>%
  mutate(
    stars = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01  ~ "**",
      P.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    label_color = ifelse(P.adj < 0.05, "white", "black")
  )

# Plot Heatmap (Lower Triangle)
p_heatmap <- ggplot(dunn_res_final, aes(x = Method1, y = Method2, fill = P.adj)) +
  geom_tile(color = "white") +
  # Custom scale to emphasize significant values (p < 0.05)
  scale_fill_gradientn(
    colors = c("darkred", "red", "orange", "white"),
    values = c(0, 0.01, 0.05, 1),
    limits = c(0, 1),
    name = "Adj. P-Value"
  ) +
  geom_text(aes(label = stars, color = label_color), size = 5, vjust = 0.7) +
  scale_color_identity() +
  scale_x_discrete(drop = FALSE) + 
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Pairwise Significance (% AUC Improvement, Dunn's Test)",
    x = NULL, y = NULL
  )

ggsave(file.path(output_plot_dir, "significance_tests.png"), plot = p_heatmap, width = 8, height = 7, dpi = 300)
cat("Significance heatmap saved (Lower Triangle).\n")


# -------------------------------------------------------------------------
# (7) Species Traits Plots (Mixed-Effects Models on RAW AUC)
# -------------------------------------------------------------------------
cat("\n###############################################\n")
cat("GENERATING TRAITS ANALYSIS PLOTS (RAW AUC, Mixed-Effects)\n")
cat("###############################################\n")

# 1. Summarize Raw AUC by Species & Method
trait_df <- final_data_auc %>%
  group_by(species, method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE), .groups = "drop") %>%
  left_join(species_map, by = c("species" = "Species"))

# --- RECODE TRAIT LABELS ---
trait_df <- trait_df %>%
  mutate(
    Prevalence.Level = recode(Prevalence.Level, "l"="Low", "m"="Medium", "h"="High"),
    Habitat = recode(Habitat, "e"="Seral", "f"="Forest"),
    Generalist.Specialist = recode(Generalist.Specialist, "g"="Generalist", "s"="Specialist"),
    Home.Range = recode(Home.Range, "s"="Small", "m"="Medium", "l"="Large")
  )

# Set Factor Levels and Order
trait_df$Prevalence.Level <- factor(trait_df$Prevalence.Level, levels=c("Low", "Medium", "High"))
trait_df$Home.Range <- factor(trait_df$Home.Range, levels=c("Small", "Medium", "Large"))
trait_df$Habitat <- factor(trait_df$Habitat, levels=c("Forest", "Seral"))
trait_df$Generalist.Specialist <- factor(trait_df$Generalist.Specialist, levels=c("Generalist", "Specialist"))

# Apply User's Fixed Algorithm Order (INCLUDING lat-long)
alg_order <- c("lat-long", "1-kmSq", "DBSC", "1to10", "2to10-sameObs", 
               "2to10", "rounded-4", "BayesOptClustGeo", "best-clustGeo")

trait_df <- trait_df %>% filter(method %in% alg_order)
trait_df$method <- factor(trait_df$method, levels = alg_order)

# Set global theme for sjPlot
sjPlot::set_theme(base = theme_classic(), axis.angle.x = 45)

# --- RUN LMER MODELS (Response = mean_auc) ---
# Using interaction (method:Trait) to get estimates for each combo
m_prev <- lmer(mean_auc ~ method:Prevalence.Level + (1|species), data=trait_df, control=lmerControl(check.rankX="silent.drop.cols"))
m_hab  <- lmer(mean_auc ~ method:Habitat + (1|species), data=trait_df, control=lmerControl(check.rankX="silent.drop.cols"))
m_spec <- lmer(mean_auc ~ method:Generalist.Specialist + (1|species), data=trait_df, control=lmerControl(check.rankX="silent.drop.cols"))
m_home <- lmer(mean_auc ~ method:Home.Range + (1|species), data=trait_df, control=lmerControl(check.rankX="silent.drop.cols"))

# Helper for plotting mixed effects with sjPlot (Plots & Lines)
plot_lmer_effects <- function(model, trait_col, title_str) {
  
  # Ensure colors match the exact factor order used in the model
  model_colors <- unname(colors[levels(trait_df$method)])
  
  sjPlot::plot_model(
    model, 
    type = "pred", 
    terms = c(trait_col, "method"), # X-axis = Trait, Group = Method
    title = title_str, 
    axis.title = c(title_str, "AUC"),
    colors = model_colors, 
    dodge = 0.6,
    legend.title = "Algorithm"
  ) 
}

# Create 4 panels
p_prev <- plot_lmer_effects(m_prev, "Prevalence.Level", "Prevalence")
p_hab  <- plot_lmer_effects(m_hab, "Habitat", "Habitat")
p_spec <- plot_lmer_effects(m_spec, "Generalist.Specialist", "Generalist/Specialist")
p_home <- plot_lmer_effects(m_home, "Home.Range", "Home Range")

# Combine with Common Legend
p_traits <- ( (p_prev + p_hab) / (p_spec + p_home) ) + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Algorithm Performance by Species Traits (Raw AUC)") & 
  theme(legend.position = "bottom")

ggsave(file.path(output_plot_dir, "traits.png"), plot = p_traits, width = 12, height = 10, dpi = 300)
cat("Traits plot saved (Mixed-Effects/sjPlot/Raw AUC).\n")


# -------------------------------------------------------------------------
# (8) Kappa Difference Plot & Saving BayesOpt Params
# -------------------------------------------------------------------------
cat("\n###############################################\n")
cat("GENERATING KAPPA DIFFERENCE PLOT\n")
cat("###############################################\n")

# 1. Load best-clustGeo params (Grid Search Best)
grid_kappa <- best_params_combined %>%
  filter(metric == "AUC") %>%
  mutate(kappa_grid = as.numeric(kappa)) %>%
  select(species, kappa_grid)

# 2. Load BayesOpt params
bayes_file <- file.path(output_dir, "BayesOptClustGeo_params.csv")
if(file.exists(bayes_file)) {
  bayes_hist <- read.csv(bayes_file)
  
  # Select best kappa per species (Maximize AUC)
  bayes_kappa <- bayes_hist %>%
    group_by(species) %>%
    slice_max(AUC, n = 1, with_ties = FALSE) %>%
    mutate(kappa_bayes = as.numeric(kappa)) %>%
    select(species, kappa_bayes, rho, AUC)
  
  # --- SAVE BAYES PARAMS SUMMARY (Requested) ---
  write.csv(bayes_kappa, file.path(output_plot_dir, "BayesOptClustGeo_best_params.csv"), row.names = FALSE)
  cat("Saved BayesOptClustGeo_best_params.csv\n")
  
  # 3. Join and Compare
  bayes_kappa_join <- bayes_kappa %>%
    left_join(species_map %>% select(Species, Abbreviation), by = c("species" = "Abbreviation")) %>%
    mutate(species_full = ifelse(!is.na(Species), Species, species)) %>%
    ungroup() %>%
    select(species_full, kappa_bayes)
  
  kappa_comp <- grid_kappa %>%
    inner_join(bayes_kappa_join, by = c("species" = "species_full")) %>%
    mutate(kappa_diff = kappa_grid - kappa_bayes)
  
  # --- SAVE COMPARISON TABLE ---
  write.csv(kappa_comp, file.path(output_plot_dir, "clustGeo_parameter_comparison.csv"), row.names = FALSE)
  
  # 4. Plot
  p_kappa <- ggplot(kappa_comp, aes(x = reorder(species, kappa_diff), y = kappa_diff)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.7) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_classic() +
    labs(
      title = "Difference in Selected Kappa",
      subtitle = NULL, 
      y = "Kappa Difference (best-clustGeo - BayesOptClustGeo)",
      x = "Species"
    )
  
  ggsave(file.path(output_plot_dir, "kappa_diff.png"), plot = p_kappa, width = 8, height = 8, dpi = 300)
  cat("Kappa difference plot saved.\n")
  
} else {
  cat("Skipping Kappa Plot: BayesOptClustGeo_params.csv not found.\n")
}


# -------------------------------------------------------------------------
# (9) Generate Maps (Psi Only) - WGS84
# -------------------------------------------------------------------------
cat("\n###############################################\n")
cat("GENERATING SPECIES MAPS (PSI ONLY - WGS84)\n")
cat("###############################################\n")

# Load Parameters
params_df <- read.csv(file.path(output_dir, "estimated_parameters.csv"))

# Load Best clustGeo Definition
if(file.exists(file.path(output_plot_dir, "best-clustGeo_params.csv"))) {
  best_cg_params <- read.csv(file.path(output_plot_dir, "best-clustGeo_params.csv"))
} else {
  best_cg_params <- NULL
}

species_list <- as.character(species_map$Abbreviation)

# Define Spatial Constants
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
wgs84_crs_str  <- "EPSG:4326"
boundary_shapefile_path <- file.path("state_covariate_raster", "boundary", "boundary.shp")

# Prepare Rasters
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
names(state_cov_raster_raw) <- c("elevation", "TCB", "TCG", "TCW", "TCA") 
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = 100)

valid_boundary <- terra::vect(boundary_shapefile_path)
valid_boundary_proj <- terra::project(valid_boundary, albers_crs_str)
cov_tif_albers_raw <- terra::mask(cov_tif_albers_raw, valid_boundary_proj)
cov_tif_albers_raw <- terra::crop(cov_tif_albers_raw, valid_boundary_proj)

# Standardize
cov_tif_albers <- cov_tif_albers_raw
for(nm in names(cov_tif_albers)) {
  mu <- mean(values(cov_tif_albers[[nm]]), na.rm=TRUE)
  sd_val <- sd(values(cov_tif_albers[[nm]]), na.rm=TRUE)
  cov_tif_albers[[nm]] <- (cov_tif_albers[[nm]] - mu) / sd_val
}
cell_area_km2 <- (100 / 1000) * (100 / 1000)

# Background for plotting
bg_raster_wgs84 <- terra::project(state_cov_raster_raw[["elevation"]], wgs84_crs_str)
valid_boundary_wgs84 <- terra::project(valid_boundary, wgs84_crs_str)
bg_raster_wgs84 <- terra::mask(bg_raster_wgs84, valid_boundary_wgs84)
bg_raster_wgs84 <- terra::crop(bg_raster_wgs84, valid_boundary_wgs84)
bbox_full <- terra::ext(bg_raster_wgs84)
bg_df_wgs84 <- as.data.frame(bg_raster_wgs84, xy = TRUE, na.rm = TRUE)

methods_for_maps <- c("1to10", "2to10", "2to10-sameObs", "1-kmSq", "lat-long", "rounded-4", "DBSC", "BayesOptClustGeo", "best-clustGeo")

predict_occuN_psi <- function(cov_stack, param_row, state_covs, cell_area) {
  intercept <- param_row$state_intercept
  betas <- as.numeric(param_row[state_covs])
  lin_pred <- cov_stack[[1]] * 0 + intercept
  for(i in seq_along(state_covs)) {
    lin_pred <- lin_pred + (cov_stack[[state_covs[i]]] * betas[i])
  }
  lambda_rast <- exp(lin_pred) * cell_area
  psi_rast <- 1 - exp(-lambda_rast)
  names(psi_rast) <- "psi"
  return(psi_rast)
}

for (sp in species_list) {
  
  cat(sprintf("Generating map for %s (WGS84)...\n", sp))
  train_file <- file.path("checklist_data", "species", sp, paste0(sp, "_zf_filtered_region_2017.csv"))
  test_file  <- file.path("checklist_data", "species", sp, paste0(sp, "_zf_filtered_region_2018.csv"))
  
  if(!file.exists(train_file) || !file.exists(test_file)) next
  
  obs_train <- read.delim(train_file, sep=",")
  obs_test  <- read.delim(test_file, sep=",")
  
  # --- FIX: Ensure compatible types before binding ---
  obs_train$observation_count <- as.character(obs_train$observation_count)
  obs_test$observation_count <- as.character(obs_test$observation_count)
  # ---------------------------------------------------
  
  obs_train <- obs_train[!is.na(obs_train$duration_minutes) & obs_train$observation_date >= "2017-05-15" & obs_train$observation_date <= "2017-07-09",]
  obs_test  <- obs_test[!is.na(obs_test$duration_minutes) & obs_test$observation_date >= "2018-05-15" & obs_test$observation_date <= "2018-07-09",]
  
  pts_df <- bind_rows(obs_train, obs_test) %>%
    mutate(
      species_observed_label = ifelse(species_observed == 1 | species_observed == TRUE, "Detection", "Non-detection"),
      species_observed_label = factor(species_observed_label, levels = c("Non-detection", "Detection"))
    ) %>% arrange(species_observed_label)
  
  obs_plot <- ggplot() + 
    geom_rect(aes(xmin = bbox_full$xmin, xmax = bbox_full$xmax, ymin = bbox_full$ymin, ymax = bbox_full$ymax), fill = "darkgray", show.legend = FALSE) +
    geom_raster(data = bg_df_wgs84, aes(x = x, y = y), fill = "#E6E6E6", show.legend = FALSE) +
    geom_point(data = pts_df, aes(x = longitude, y = latitude, color = species_observed_label, shape = species_observed_label, fill = species_observed_label, size = species_observed_label), show.legend = TRUE) +
    scale_color_manual(name = "Observation", values = c("Detection" = "black", "Non-detection" = "black")) +
    scale_fill_manual(name = "Observation", values = c("Detection" = "#39FF14", "Non-detection" = "#83A1CD")) +
    scale_shape_manual(name = "Observation", values = c("Detection" = 24, "Non-detection" = 22)) +
    scale_size_manual(name = "Observation", values = c("Detection" = 1.6, "Non-detection" = 1.5)) +
    labs(title = "Species Observations") +
    theme_void() +
    coord_fixed(ratio = 1.0, xlim = c(bbox_full$xmin, bbox_full$xmax), ylim = c(bbox_full$ymin, bbox_full$ymax), expand = FALSE) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = -25, face = "bold", size = 18),
      legend.position = "inside", legend.position.inside = c(0.5, -0.125), legend.direction = "vertical",
      legend.text = element_text(size = 14), legend.title = element_text(size = 14, margin = margin(l = 20, b = 10)), legend.spacing.x = unit(3, "cm"), legend.key.size = unit(1, 'cm')
    )
  
  psi_plots <- list()
  sp_params <- params_df %>% filter(species == sp)
  
  full_name_row <- species_map %>% filter(Abbreviation == sp)
  full_name <- if(nrow(full_name_row) > 0) full_name_row$Species[1] else sp
  
  best_cg_row <- best_cg_params %>% filter(species == full_name, metric == "AUC")
  actual_best_method <- if(nrow(best_cg_row) > 0) best_cg_row$method[1] else NA
  
  for (i in seq_along(methods_for_maps)) {
    m_label <- methods_for_maps[i]
    m_lookup <- if (m_label == "best-clustGeo") actual_best_method else m_label
    plot_title <- if (m_label == "best-clustGeo") "best-clustGeo" else m_label
    
    m_param <- sp_params %>% filter(method == m_lookup)
    
    if (nrow(m_param) == 0 || is.na(m_lookup)) {
      psi_plots[[i]] <- ggplot() + theme_void() + labs(title = plot_title) + theme(plot.title = element_text(hjust = 0.5, size = 14))
      next
    }
    
    psi_rast_albers <- predict_occuN_psi(cov_tif_albers, m_param, c("elevation","TCB","TCG","TCW","TCA"), cell_area_km2)
    psi_rast_wgs84 <- terra::project(psi_rast_albers, wgs84_crs_str)
    psi_rast_wgs84 <- terra::mask(psi_rast_wgs84, valid_boundary_wgs84)
    psi_rast_wgs84 <- terra::crop(psi_rast_wgs84, valid_boundary_wgs84)
    psi_df <- as.data.frame(psi_rast_wgs84, xy = TRUE, na.rm=TRUE)
    
    psi_plots[[i]] <- ggplot() +
      geom_rect(aes(xmin = bbox_full$xmin, xmax = bbox_full$xmax, ymin = bbox_full$ymin, ymax = bbox_full$ymax), fill = "darkgray", show.legend = FALSE) +
      geom_raster(data = psi_df, aes(x = x, y = y, fill = psi)) +
      scale_fill_viridis_c(option = "B", limits = c(0.0, 1.0), name = "Occupancy Probability") +
      theme_void() + 
      coord_fixed(ratio = 1.0, xlim = c(bbox_full$xmin, bbox_full$xmax), ylim = c(bbox_full$ymin, bbox_full$ymax), expand = FALSE) +
      labs(title = plot_title) +
      theme(
        legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 14, vjust = 1), 
        legend.key.width = unit(1, "cm"), legend.box.margin = margin(t = 20), plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
  }
  
  grid_p <- ggarrange(plotlist = psi_plots, nrow = 2, ncol = 5, common.legend = TRUE, legend = "bottom")
  final <- (obs_plot + grid_p + plot_layout(nrow = 1, widths = c(1, 3.5)))
  ggsave(file.path(map_output_dir, paste0(sp, ".png")), plot = final, width = 17, height = 9.5, dpi = 300)
}

cat("Maps generated.\n")