# -----------------------------------------------------------------
# Real Species Experiments for occuN model
# -----------------------------------------------------------------

###
# 1. SETUP
###

install_now = FALSE
if (install_now){
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  if (!requireNamespace("devtools", quietly = FALSE)) install.packages("devtools")
  suppressMessages(devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE))
}

library(unmarked)
library(dplyr)
library(tidyr)
library(PRROC)
library(terra)

source(file.path("R", "utils.R"))
source(file.path("R", "data_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))
source(file.path("R", "plotting_helpers.R"))

set.seed(123) 

###
# 2. CONFIGS
###

# Species list
species_names <- c(
    "AMCR", "AMRO", "BAEA", "BKHGRO", "BRCR", "BUTI", "CASC", "CHBCHI", 
    "COHA", "HAFL", "HAWO", "HEWA", "MAWA", "MOQU", "NOFL", "NOOW", 
    "OLFL", "PAFL", "PAWR", "PIWO", "REHA", "SOSP", "SPTO", "SWTH", 
    "WAVI", "WEPE", "WETA", "WIWA", "WRENTI", "YEBCHA", "YEWA"
)

# Comparison methods


method_names <- c(
    "1to10", 
    "2to10", 
    "2to10-sameObs", 
    "1-kmSq",
    "lat-long", 
    "rounded-4", 
    "SVS", 
    "1-per-UL",
    "clustGeo-25-60",
    "clustGeo-50-60",
    "clustGeo-75-60",
    "clustGeo-25-70",
    "clustGeo-50-70",
    "clustGeo-75-70",
    "clustGeo-25-80",
    "clustGeo-50-80",
    "clustGeo-75-80",
    "clustGeo-25-90",
    "clustGeo-50-90",
    "clustGeo-75-90",
    "DBSC",
    "BayesOptClustGeo"
)

# Methods to plot
methods_to_plot <- c(
    "1to10", "2to10", "2to10-sameObs", "lat-long", "SVS", "1-per-UL",
    "1-kmSq", "rounded-4", "DBSC", "BayesOptClustGeo"
)

methods_to_plot_clustGeo <- c(
    "clustGeo-25-60", "clustGeo-50-60", "clustGeo-75-60", "clustGeo-25-70", "clustGeo-50-70", "clustGeo-75-70",
    "clustGeo-25-80", "clustGeo-50-80", "clustGeo-75-80", "clustGeo-25-90", "clustGeo-50-90", "clustGeo-75-90"
)


# Covariates
state_cov_names <- c("elevation", "TCB", "TCG", "TCW", "TCA")
obs_cov_names <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")

# Optimization & Simulation Settings
selected_optimizer <- "nlminb"
n_fit_repeats <- 100
n_test_repeats <- 25

res_m <- 100 
buffer_m <- 200

PARAM_LOWER <- -10
PARAM_UPPER <- 10
INIT_LOWER <- -2
INIT_UPPER <- 2

# Output Directory
output_dir <- file.path("species_experiments", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

###
# 3. PREPROCESS RASTER (Global)
###

cat("--- Pre-processing global raster data... ---\n")

# 1. Define CRS
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# 2. Load Native Raster
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster_raw) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster_raw) <- state_cov_names

# 3. Project to Albers (Raw Values)
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = res_m)

# 4. STANDARDIZE RASTER
standardization_results <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- standardization_results$raster
state_cov_params <- standardization_results$params

# 5. Area Raster
cell_area_km2 <- (res_m / 1000) * (res_m / 1000)
area_j_raster <- cov_tif_albers[[1]] * 0 + cell_area_km2
names(area_j_raster) <- "area"

# 6. Full Raster Covs (for occuN preparation)
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))[, state_cov_names, drop = FALSE]
full_raster_covs[is.na(full_raster_covs)] <- 0

# Boundary for plotting
boundary_shapefile_path <- file.path("state_covariate_raster", "boundary", "boundary.shp")


###
# 4. MAIN LOOP (Per Species)
###

all_param_results <- list()
all_pred_results <- list()
all_clustering_stats_pre <- list()
all_clustering_stats_post <- list()

for (species_name in species_names) {
  
  cat(sprintf("\n\n###############################################\n"))
  cat(sprintf("PROCESSING SPECIES: %s\n", species_name))
  cat(sprintf("###############################################\n"))
  
  # === 4.1 DATA PREPARATION ===
  cat("--- Preparing Train Data ---\n")
  
  # New prepare_train_data preserves species_observed.
  train_data_res <- prepare_train_data(
    state_covs = state_cov_names, 
    obs_covs = obs_cov_names, 
    cov_tif = cov_tif_albers_raw, 
    state_standardization_params = state_cov_params,
    placeholder_spec_name = species_name
  )
  
  base_train_df <- train_data_res$train_df
  full_standardization_params <- train_data_res$standardization_params
  
  cat("--- Preparing Test Data ---\n")
  base_test_df <- prepare_test_data(
    state_covs = state_cov_names, 
    obs_covs = obs_cov_names, 
    cov_tif = cov_tif_albers_raw, 
    standardization_params = full_standardization_params,
    placeholder_spec_name = species_name
  )

  # === 4.2 CLUSTERING (TRAIN) ===
  cat("--- Computing clusterings... ---\n")
  all_clusterings <- get_clusterings(method_names, base_train_df, state_cov_names, NULL)
  
  cat("--- Computing initial site geometries... ---\n")
  all_site_geometries <- list()
  for (method_name in method_names) {
    cluster_data <- all_clusterings[[method_name]]
    if (is.list(cluster_data) && "result_df" %in% names(cluster_data)) cluster_data <- cluster_data$result_df
    
    if (!is.null(cluster_data)) {
      all_site_geometries[[method_name]] <- create_site_geometries(cluster_data, cov_tif_albers, buffer_m, method_name)
    }
  }
  
  # Record Stats PRE Split
  stats_pre <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
  stats_pre$species <- species_name
  all_clustering_stats_pre[[length(all_clustering_stats_pre) + 1]] <- stats_pre

  # --- PLOTTING PRE-SPLIT ---
  cat("--- Plotting sites (PRE-SPLIT)... ---\n")
  try({
    plot_sites(
      base_train_df = base_train_df,
      all_clusterings = all_clusterings,
      all_site_geometries = all_site_geometries,
      elevation_raster = cov_tif_albers_raw, 
      methods_to_plot = methods_to_plot,
      boundary_shp_path = boundary_shapefile_path,
      output_path = file.path(output_dir, paste0(species_name, "_site_visualization_PRE.png")),
      cluster_labels = TRUE
    )
  })
  
  
  # === 4.3 SPLIT DISJOINT SITES ===
  cat("--- Post-processing: Splitting disjoint geometries... ---\n")
  for (method_name in names(all_site_geometries)) {
    curr_geoms <- all_site_geometries[[method_name]]
    curr_data_obj <- all_clusterings[[method_name]]
    
    is_list_obj <- is.list(curr_data_obj) && "result_df" %in% names(curr_data_obj)
    curr_data <- if(is_list_obj) curr_data_obj$result_df else curr_data_obj
    
    split_res <- disjoint_site_geometries(curr_geoms, curr_data)
    all_site_geometries[[method_name]] <- split_res$geoms
    
    if (is_list_obj) {
      all_clusterings[[method_name]]$result_df <- split_res$data
    } else {
      all_clusterings[[method_name]] <- split_res$data
    }
  }
  
  # Record Stats POST Split
  stats_post <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
  stats_post$species <- species_name
  all_clustering_stats_post[[length(all_clustering_stats_post) + 1]] <- stats_post
  
  # --- PLOTTING POST-SPLIT ---
  cat("--- Plotting sites (POST-SPLIT)... ---\n")
  try({
    plot_sites(
      base_train_df = base_train_df,
      all_clusterings = all_clusterings,
      all_site_geometries = all_site_geometries,
      elevation_raster = cov_tif_albers_raw, 
      methods_to_plot = methods_to_plot,
      boundary_shp_path = boundary_shapefile_path,
      output_path = file.path(output_dir, paste0(species_name, "_site_visualization_POST.png")),
      cluster_labels = TRUE
    )
  })
  
  # === 4.4 W MATRICES ===
  cat("--- Generating W matrices... ---\n")
  all_w_matrices <- list()
  for (m_name in names(all_site_geometries)) {
    if (!is.null(all_site_geometries[[m_name]])) {
      all_w_matrices[[m_name]] <- generate_overlap_matrix(all_site_geometries[[m_name]], cov_tif_albers)
    }
  }
  
  
  # === 4.5 TEST SPATIAL STRUCTURES ===
  cat("--- Preparing Test Structures ---\n")
  test_structures <- prepare_test_spatial_structures(
    test_df = base_test_df,
    albers_crs = albers_crs_str,
    buffer_m = buffer_m,
    cov_raster_albers = cov_tif_albers,
    area_raster = area_j_raster
  )
  base_test_df_ready <- test_structures$test_df
  w_matrix_test <- test_structures$w_matrix
  
  # Pre-generate Test Splits
  test_splits_list <- list()
  for (r in 1:n_test_repeats) {
    test_splits_list[[r]] <- spatial_subsample_dataset(base_test_df_ready, res_m/1000, r)
  }

  
  # === 4.6 MODEL FITTING LOOP ===
  for (method_name in method_names) {
    cat(sprintf("\n  [Species: %s] Fitting Method: %s\n", species_name, method_name))
    
    # Prep occuN Data
    current_clustering_df <- all_clusterings[[method_name]]
    if (is.list(current_clustering_df) && "result_df" %in% names(current_clustering_df)) {
      current_clustering_df <- current_clustering_df$result_df
    }
    
    w_matrix <- all_w_matrices[[method_name]]
    
    if (is.null(w_matrix)) {
       cat(sprintf("    Skipping %s (No W matrix)\n", method_name)); next
    }
    
    # Use modular function to prep UMF
    # Note: Dummy site column is no longer needed.
    umf <- prepare_occuN_data(base_train_df, current_clustering_df, w_matrix, obs_cov_names, full_raster_covs)
    
    obs_formula <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
    state_formula <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))
    
    cat(sprintf("    Fitting %s (M=%d)... ", method_name, nrow(umf@y)))
    
    # Fit Model
    fm <- fit_occuN_model(
      umf, state_formula, obs_formula,
      n_reps = n_fit_repeats, stable_reps = n_fit_repeats,
      optimizer = selected_optimizer, lower = PARAM_LOWER, upper = PARAM_UPPER,
      init_lower = INIT_LOWER, init_upper = INIT_LOWER
    )
    
    if (is.null(fm)) {
      cat("FAILED.\n")
      # Record Failure (optional, or just skip)
      next
    }
    cat("Done.\n")
    
    # --- SAVE PARAMETERS ---
    est_alphas <- coef(fm, 'det')   
    est_betas <- coef(fm, 'state')  
    
    state_col_names <- c("state_intercept", state_cov_names)
    obs_col_names   <- c("obs_intercept", obs_cov_names)
    
    param_row <- data.frame(
      species = species_name,
      method = method_name,
      nll = fm@negLogLike,
      convergence = 0
    )
    
    for(i in seq_along(state_col_names)) param_row[[state_col_names[i]]] <- est_betas[i]
    for(i in seq_along(obs_col_names)) param_row[[obs_col_names[i]]] <- est_alphas[i]
    
    all_param_results[[length(all_param_results) + 1]] <- param_row
    
    
    # --- PREDICT & TEST (Repeats) ---
    for (repeat_num in 1:n_test_repeats) {
      test_df <- test_splits_list[[repeat_num]]
      
      X_state <- model.matrix(state_formula, data = test_df)
      X_obs <- model.matrix(obs_formula, data = test_df)
      
      # Predictions
      pred_psi <- 1 - exp(-(exp(X_state %*% est_betas) * test_df$area_j))
      pred_det <- plogis(X_obs %*% est_alphas)
      pred_obs_prob <- pred_psi * pred_det 
      
      metrics <- calculate_classification_metrics(pred_obs_prob, test_df$species_observed)
      
      pred_row <- data.frame(
        species = species_name,
        method = method_name,
        test_repeat = repeat_num,
        auc = metrics$auc,
        auprc = metrics$auprc
      )
      all_pred_results[[length(all_pred_results) + 1]] <- pred_row
    }
    
    # Memory Cleanup (per method)
    rm(umf, fm)
    gc()
    
  } # End Method Loop
  
  # Memory Cleanup (per species)
  rm(all_clusterings, all_site_geometries, all_w_matrices, test_structures, test_splits_list)
  gc()
  
} # End Species Loop

# Save Final Results
write.csv(dplyr::bind_rows(all_param_results), file.path(output_dir, "estimated_parameters.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_pred_results), file.path(output_dir, "predictive_performance.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_clustering_stats_pre), file.path(output_dir, "clustering_stats_PRE.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_clustering_stats_post), file.path(output_dir, "clustering_stats_POST.csv"), row.names = FALSE)

cat("Done.\n")