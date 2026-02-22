#######################################
# Plot Species Cluster-Based Experiments

# February 20, 2026
#######################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# Create output directory
# dir.create("simulation_experiments/output/plots", recursive = TRUE, showWarnings = FALSE)

output_dir <- file.path("output", "species_experiments", "clusters")
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
# We exclude clustGeo variants and SLIC variants from this set
other_data <- data %>% 
  filter((!grepl("^clustGeo", method) | method == "BayesOptClustGeo"))

# --- Process ClustGeo ---
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


# --- Combine all data ---
final_data_auc <- bind_rows(other_data, best_cg_data_auc)
final_data_auprc <- bind_rows(other_data, best_cg_data_auprc)


# -------------------------------------------------------------------------
# (1.5) Save Best clustGeo Parameters to CSV
# -------------------------------------------------------------------------

# Combine the best methods identified for AUC and AUPRC
best_params_auc <- best_cg_auc_methods %>%
  mutate(metric = "AUC") %>%
  rename(method = best_method_auc)

best_params_auprc <- best_cg_auprc_methods %>%
  mutate(metric = "AUPRC") %>%
  rename(method = best_method_auprc)

best_params_combined <- bind_rows(best_params_auc, best_params_auprc) %>%
  # Parse the method string "clustGeo-rho-kappa"
  # Removes "clustGeo-" and splits the remaining numbers
  mutate(params_str = sub("^clustGeo-", "", method)) %>%
  separate(params_str, c("rho", "kappa"), sep = "-") %>%
  # Select and reorder final columns
  select(species, metric, method, rho, kappa)

# Save to the plots directory
write.csv(best_params_combined, file.path(output_plot_dir, "best-clustGeo_params.csv"), row.names = FALSE)

cat("Best clustGeo parameters saved to:", file.path(output_plot_dir, "best-clustGeo_params.csv"), "\n")

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
               "rounded-4", "best-clustGeo", "DBSC", "BayesOptClustGeo")
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



# -------------------------------------------------------------------------
# (5) Generate Maps (Psi Only) - WGS84 (Final Layout)
# -------------------------------------------------------------------------

cat("\n###############################################\n")
cat("GENERATING SPECIES MAPS (PSI ONLY - WGS84)\n")
cat("###############################################\n")

# --- 0. SETUP & LIBRARIES ---
library(patchwork)  # For layout
library(terra)      # For raster operations
library(sf)         # For vector operations
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Source helpers to access 'standardize_state_covs' if available
if(file.exists(file.path("R", "data_helpers.R"))) source(file.path("R", "data_helpers.R"))

# 1. Setup Output Directory
map_output_dir <- file.path(output_plot_dir, "maps")
if (!dir.exists(map_output_dir)) dir.create(map_output_dir, recursive = TRUE)

# 2. Load Parameters
params_df <- read.csv(file.path(output_dir, "estimated_parameters.csv"))

# --- Load Best clustGeo Definition ---
best_cg_file <- file.path(output_plot_dir, "best-clustGeo_params.csv")
if(file.exists(best_cg_file)) {
  best_cg_params <- read.csv(best_cg_file)
} else {
  warning("best-clustGeo_params.csv not found. 'best-clustGeo' map may fail.")
  best_cg_params <- NULL
}

# --- Define species_names from the already loaded species_map ---
species_names <- as.character(species_map$Abbreviation)

# 3. Define Spatial Constants & Load Raster
# Model uses Albers, Plotting uses WGS84
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
wgs84_crs_str  <- "EPSG:4326"
boundary_shapefile_path <- file.path("state_covariate_raster", "boundary", "boundary.shp")

cat("Loading and preparing environmental rasters...\n")

# A. Prepare ALBERS Raster (For Model Predictions)
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
names(state_cov_raster_raw) <- c("elevation", "TCB", "TCG", "TCW", "TCA") 
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = 100)

# Mask and Crop to Boundary (in Albers)
valid_boundary <- terra::vect(boundary_shapefile_path)
valid_boundary_proj <- terra::project(valid_boundary, albers_crs_str)
cov_tif_albers_raw <- terra::mask(cov_tif_albers_raw, valid_boundary_proj)
cov_tif_albers_raw <- terra::crop(cov_tif_albers_raw, valid_boundary_proj)

# Standardize Albers Raster
if(exists("standardize_state_covs")) {
  standardization_results <- standardize_state_covs(cov_tif_albers_raw)
  cov_tif_albers <- standardization_results$raster
} else {
  cat("Helper not found, performing manual standardization...\n")
  cov_tif_albers <- cov_tif_albers_raw
  for(nm in names(cov_tif_albers)) {
    mu <- mean(values(cov_tif_albers[[nm]]), na.rm=TRUE)
    sd_val <- sd(values(cov_tif_albers[[nm]]), na.rm=TRUE)
    cov_tif_albers[[nm]] <- (cov_tif_albers[[nm]] - mu) / sd_val
  }
}
cell_area_km2 <- (100 / 1000) * (100 / 1000)

# B. Prepare WGS84 Background (For Plotting)
bg_raster_wgs84 <- terra::project(state_cov_raster_raw[["elevation"]], wgs84_crs_str)
valid_boundary_wgs84 <- terra::project(valid_boundary, wgs84_crs_str)
bg_raster_wgs84 <- terra::mask(bg_raster_wgs84, valid_boundary_wgs84)
bg_raster_wgs84 <- terra::crop(bg_raster_wgs84, valid_boundary_wgs84)

# Get Extent for coord_fixed
bbox_full <- terra::ext(bg_raster_wgs84)
bg_df_wgs84 <- as.data.frame(bg_raster_wgs84, xy = TRUE, na.rm = TRUE)

# 4. Define Methods to Map
methods_for_maps <- c(
  "1to10", 
  "2to10", 
  "2to10-sameObs", 
  "1-kmSq",
  "lat-long", 
  "rounded-4", 
  "DBSC", 
  "BayesOptClustGeo",
  "best-clustGeo" 
)

# 5. Define Prediction Function (Psi Only)
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

# 6. Loop over Species
for (sp in species_names) {
  
  cat(sprintf("Generating map for %s (WGS84)...\n", sp))
  
  # --- A. Prepare Observations ---
  train_file <- file.path("checklist_data", "species", sp, paste0(sp, "_zf_filtered_region_2017.csv"))
  test_file  <- file.path("checklist_data", "species", sp, paste0(sp, "_zf_filtered_region_2018.csv"))
  
  if(!file.exists(train_file) || !file.exists(test_file)) {
    cat(sprintf("Skipping %s (Data files not found)\n", sp))
    next
  }

  obs_train <- read.delim(train_file, sep=",")
  obs_test  <- read.delim(test_file, sep=",")
  
  if("observation_count" %in% names(obs_train)) obs_train$observation_count <- as.character(obs_train$observation_count)
  if("observation_count" %in% names(obs_test))  obs_test$observation_count  <- as.character(obs_test$observation_count)

  obs_train <- obs_train[!is.na(obs_train$duration_minutes) & obs_train$observation_date >= "2017-05-15" & obs_train$observation_date <= "2017-07-09",]
  obs_test  <- obs_test[!is.na(obs_test$duration_minutes) & obs_test$observation_date >= "2018-05-15" & obs_test$observation_date <= "2018-07-09",]
  
  pts_df <- bind_rows(obs_train, obs_test) %>%
    mutate(
      species_observed_label = ifelse(species_observed == 1 | species_observed == TRUE, "Detection", "Non-detection"),
      species_observed_label = factor(species_observed_label, levels = c("Non-detection", "Detection"))
    ) %>%
    arrange(species_observed_label)
  
  # --- B. Create Left Plot (Observations - WGS84) ---
  
  obs_plot <- ggplot() + 
    # Background Box
    geom_rect(aes(xmin = bbox_full$xmin, xmax = bbox_full$xmax, 
                  ymin = bbox_full$ymin, ymax = bbox_full$ymax), 
              fill = "darkgray", show.legend = FALSE) +
    
    # Landmass (Solid Light Grey)
    geom_raster(data = bg_df_wgs84, aes(x = x, y = y), fill = "#E6E6E6", show.legend = FALSE) +
    
    # Observation Points
    geom_point(data = pts_df, aes(x = longitude, y = latitude, 
                                  color = species_observed_label, 
                                  shape = species_observed_label, 
                                  fill = species_observed_label, 
                                  size = species_observed_label), show.legend = TRUE) +
    
    scale_color_manual(name = "Observation", values = c("Detection" = "black", "Non-detection" = "black")) +
    scale_fill_manual(name = "Observation", values = c("Detection" = "#39FF14", "Non-detection" = "#83A1CD")) +
    scale_shape_manual(name = "Observation", values = c("Detection" = 24, "Non-detection" = 22)) +
    scale_size_manual(name = "Observation", values = c("Detection" = 1.6, "Non-detection" = 1.5)) +
    
    labs(title = "Species Observations") +
    theme_void() +
    coord_fixed(
        ratio = 1.0, 
        xlim = c(bbox_full$xmin, bbox_full$xmax), 
        ylim = c(bbox_full$ymin, bbox_full$ymax), 
        expand = FALSE
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold", size = 24), # Increased Title Size
      legend.position = "inside",
      legend.position.inside = c(0.5, -0.16),
      legend.direction = "vertical",
      legend.text = element_text(size = 14), 
      legend.title = element_text(size = 15, margin = margin(l = 20, b = 10)), 
      legend.spacing.x = unit(3, "cm"), 
      legend.key.size = unit(1, 'cm'),
      plot.margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "cm")
    )

  # --- C. Generate Prediction Rasters ---
  psi_plots <- list()
  sp_params <- params_df %>% filter(species == sp)
  
  full_name_row <- species_map %>% filter(Abbreviation == sp)
  if(nrow(full_name_row) > 0) {
      full_name <- full_name_row$Species[1]
  } else {
      full_name <- sp 
  }

  best_cg_row <- best_cg_params %>% filter(species == full_name, metric == "AUC")
  if(nrow(best_cg_row) == 0) best_cg_row <- best_cg_params %>% filter(species == full_name) 
  
  actual_best_method <- if(nrow(best_cg_row) > 0) best_cg_row$method[1] else NA
  
  for (i in seq_along(methods_for_maps)) {
    m_label <- methods_for_maps[i]
    
    if (m_label == "best-clustGeo") {
      if(is.na(actual_best_method)) {
         psi_plots[[i]] <- ggplot() + theme_void() + labs(title = "best-clustGeo (NA)") + theme(plot.title = element_text(hjust = 0.5))
         next
      }
      m_lookup <- actual_best_method
      plot_title <- "best-clustGeo"
    } else {
      m_lookup <- m_label
      plot_title <- m_label
    }

    m_param <- sp_params %>% filter(method == m_lookup)
    
    if (nrow(m_param) == 0) {
      psi_plots[[i]] <- ggplot() + theme_void() + labs(title = plot_title) + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      next
    }
    
    psi_rast_albers <- predict_occuN_psi(cov_tif_albers, m_param, c("elevation","TCB","TCG","TCW","TCA"), cell_area_km2)
    psi_rast_wgs84 <- terra::project(psi_rast_albers, wgs84_crs_str)
    psi_rast_wgs84 <- terra::mask(psi_rast_wgs84, valid_boundary_wgs84)
    psi_rast_wgs84 <- terra::crop(psi_rast_wgs84, valid_boundary_wgs84)
    psi_df <- as.data.frame(psi_rast_wgs84, xy = TRUE, na.rm=TRUE)
    
    psi_plots[[i]] <- ggplot() +
      geom_rect(aes(xmin = bbox_full$xmin, xmax = bbox_full$xmax, 
                    ymin = bbox_full$ymin, ymax = bbox_full$ymax), 
                fill = "darkgray", show.legend = FALSE) +
      geom_raster(data = psi_df, aes(x = x, y = y, fill = psi)) +
      scale_fill_viridis_c(option = "B", limits = c(0.0, 1.0), name = "Occupancy Probability") +
      theme_void() + 
      coord_fixed(
        ratio = 1.0, 
        xlim = c(bbox_full$xmin, bbox_full$xmax), 
        ylim = c(bbox_full$ymin, bbox_full$ymax), 
        expand = FALSE
      ) +
      labs(title = plot_title) +
      theme(
        legend.position = "bottom", 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16, vjust = 1), 
        legend.key.width = unit(1.5, "cm"),
        legend.box.margin = margin(t = 30), # Move legend down
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold") # Bigger Title
      )
  }
  
  # --- D. Assemble and Save ---
  grid_p <- ggarrange(plotlist = psi_plots, nrow = 2, ncol = 5, common.legend = TRUE, legend = "bottom")
  
  # Adjusted widths: 1:4 makes the left plot smaller relative to the previous 1:2.5
  final <- obs_plot + grid_p + plot_layout(nrow = 1, widths = c(1, 4))
  
  ggsave(file.path(map_output_dir, paste0(sp, ".png")), plot = final, width = 17, height = 9.5, dpi = 300)
}

cat("Maps generated successfully in output/species_experiments/clusters/plots/maps/\n")