# -----------------------------------------------------------------
# Simulation for occuN model
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
source(file.path("R", "simulation_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))
source(file.path("R", "plotting_helpers.R"))

set.seed(123) 

###
# 2. CONFIGS
###
comparison_method_list <- c(
  "1to10", "2to10", "2to10-sameObs", "1-kmSq",
  "lat-long", "rounded-4", "SVS", "1-per-UL",
  "DBSC", "BayesOptClustGeo"
)

comparison_method_list <- c(
  "1-kmSq"
)


selected_optimizer <- "nlminb"

sim_params <- read.delim(file.path("config", "simulation_parameters.csv"), sep = ",", header = T)
sim_clusterings <- read.delim(file.path("config", "simulation_clusterings.csv"), sep = ",", header = T)

###
# 3. SIMULATION SETTINGS
###
n_simulations <- 25
n_fit_repeats <- 25
n_test_repeats <- 25

# n_simulations <- 1 # Debug override
# n_fit_repeats <- 25 # Debug override
# n_test_repeats <- 1 # Debug override


res_m <- 100 
buffer_m <- 200

state_cov_names <- names(sim_params)[2:6]
obs_cov_names <- names(sim_params)[8:12]




###
# 4. PREPROCESS RASTER DATA
###


# 1. Define CRS
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# 2. Load Native Raster
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster_raw) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster_raw) <- state_cov_names

# 3. Project to Albers (Raw Values)
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = res_m)

# 4. STANDARDIZE RASTER FIRST (New Step)
#    Get the standardized raster AND the global parameters
standardization_results <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- standardization_results$raster
state_cov_params <- standardization_results$params

# 5. Prepare Training Data
#    Pass RAW raster for extraction (so values match the logic of applying params), 
#    BUT pass the 'state_cov_params' so the DF is scaled using Global stats.
base_train_data <- prepare_train_data(
    state_covs = state_cov_names, 
    obs_covs = obs_cov_names, 
    cov_tif = cov_tif_albers_raw, # Extract raw values
    state_standardization_params = state_cov_params # Apply global scaling
)

base_train_df <- base_train_data$train_df
# This list now has Raster Stats (for State) and Training Stats (for Obs)
full_standardization_params <- base_train_data$standardization_params 

# 6. Prepare Test Data
#    Pass the full parameter list so Test data is scaled exactly like Train data
base_test_df <- prepare_test_data(
    state_covs = state_cov_names, 
    obs_covs = obs_cov_names, 
    cov_tif = cov_tif_albers_raw, 
    standardization_params = full_standardization_params
)

area_j_raster <- cov_tif_albers[[1]] * 0 + 1
names(area_j_raster) <- "area"

# 7. Use the Standardized Raster for generating Simulation Truth
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))[, state_cov_names, drop = FALSE]
full_raster_covs[is.na(full_raster_covs)] <- 0


boundary_shapefile_path <- file.path("state_covariate_raster", "boundary", "boundary.shp")




###
# 5. TRAIN SITE GEOMETRIES
###
reference_method_list <- sim_clusterings$method
all_method_names <- unique(c(reference_method_list, comparison_method_list))

cat("--- Pre-computing ALL clusterings... ---\n")
# Pass 'cov_tif_albers' values/names depending on what clustGeo needs (it needs names + df values)
# The get_clusterings function uses the DF, so this is fine.
all_clusterings <- get_clusterings(all_method_names, base_train_df, state_cov_names, NULL)

cat("--- Pre-computing site geometries... ---\n")
all_site_geometries <- list() 
for (method_name in all_method_names) {
  cluster_data <- all_clusterings[[method_name]]
  if (is.list(cluster_data) && "result_df" %in% names(cluster_data)) cluster_data <- cluster_data$result_df
  
  if (!is.null(cluster_data)) {
    # --- CRITICAL FIX: Pass the ALBERS raster to ensure grid alignment ---
    all_site_geometries[[method_name]] <- create_site_geometries(cluster_data, cov_tif_albers, buffer_m, method_name, "km")
  }
}



mean(base_train_df$elevation)

###
# 6. CLUSTERING DESCRIPTIVE STATS
###
clustering_summary_df <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.csv(clustering_summary_df, file.path(output_dir, "clustering_descriptive_stats.csv"), row.names = FALSE)


###
# 7. CLUSTERING SIMILARITY METRICS
###

cat("--- Calculating pairwise clustering similarity metrics (ARI, AMI, NID) ---\n")

clustering_similarity_list <- list()

ref_methods_to_check <- sim_clusterings$method
comp_methods_to_check <- unique(c(sim_clusterings$method, comparison_method_list))

clustering_similarity_df <- compute_clustering_similarity( all_clusterings = all_clusterings, ref_methods = ref_methods_to_check, comp_methods = comp_methods_to_check)

write.csv(clustering_similarity_df, file.path(output_dir, "clustering_similarity_stats.csv"), row.names = FALSE)
cat(sprintf("--- Clustering similarity stats saved to %s/clustering_similarity_stats.csv ---\n", output_dir))




###
# 8. PLOT SITES
###
# all_method_names_plot_order <- c(
#   "1to10", "2to10", "2to10-sameObs", "lat-long", "SVS", "1-per-UL",
#   "0.125-kmSq", "0.25-kmSq", "0.5-kmSq", "1-kmSq", "2-kmSq", "rounded-4",
#   "clustGeo-50-20", "clustGeo-50-40", "clustGeo-50-60", "clustGeo-50-80",  
#   "BayesOptClustGeo", "DBSC"
# )

all_method_names_plot_order <- c(
  "1to10", "2to10", "2to10-sameObs", "lat-long", "SVS", "1-per-UL",
  "0.125-kmSq", "1-kmSq", "clustGeo-50-60", "BayesOptClustGeo", "DBSC", "rounded-4"
)

# Plot using the UNSCALED 100m raster
site_plot <- plot_sites(
  base_train_df = base_train_df,
  all_clusterings = all_clusterings,
  all_site_geometries = all_site_geometries,
  elevation_raster = cov_tif_albers_raw, 
  methods_to_plot = all_method_names_plot_order,
  boundary_shp_path = boundary_shapefile_path,
  output_path = file.path(output_dir, "site_cluster_visualization.png")
)


###
# 9. MEMORY CLEANUP
###
cat("--- Cleaning up heavy raster objects ---\n")

# Remove the raw unscaled raster to free memory
rm(state_cov_raster_raw, cov_tif_albers_raw)

# Remove the heavy geometry objects (we only need the W matrices now)
# (Extract W matrices first if you haven't already done so in previous steps)
all_w_matrices <- list()
for (m_name in names(all_site_geometries)) {
  if (!is.null(all_site_geometries[[m_name]])) {
    all_w_matrices[[m_name]] <- attr(all_site_geometries[[m_name]], "w_matrix")
  }
}
rm(all_site_geometries)

# Force garbage collection
gc()


###
# 10. TEST SITE GEOMETRIES
###
test_structures <- prepare_test_spatial_structures(
  test_df = base_test_df,
  albers_crs = albers_crs_str,
  buffer_m = buffer_m,
  cov_raster_albers = cov_tif_albers,
  area_raster = area_j_raster
)

base_test_df <- test_structures$test_df
w_matrix_test <- test_structures$w_matrix

rm(test_structures)
gc()


###
# 11. MAIN SIMULATION LOOP
###
all_dataset_stats <- list()
all_param_results <- list()
all_pred_results <- list()



for (cluster_idx in seq_len(nrow(sim_clusterings))) {
  current_clustering_method <- sim_clusterings$method[cluster_idx]
  cat(paste("\n=== REF CLUSTERING:", current_clustering_method, "===\n"))
  
  current_reference_dataframe <- all_clusterings[[current_clustering_method]]

  for (param_idx in seq_len(nrow(sim_params))) {
    current_parameter_set <- sim_params[param_idx, ]
    cat(paste("  Param Set:", param_idx, "\n"))

    # --- 6.1 Pre-calc N_j_raster ---
    state_par_list <- as.list(current_parameter_set[, c("state_intercept", state_cov_names)])
    names(state_par_list)[1] <- "intercept"
    
    # Calculate Lambda (Density) Vector
    log_lambda_j <- cov_tif_albers[[1]] * 0 + state_par_list$intercept
    for (nm in state_cov_names) log_lambda_j <- log_lambda_j + (cov_tif_albers[[nm]] * state_par_list[[nm]])
    
    N_j_raster <- exp(log_lambda_j) * area_j_raster
    # DENSITY VECTOR (for Matrix Mult)
    # Note: Do not multiply by area here; W matrix has area.
    # cell_density_vector <- terra::values(exp(log_lambda_j))
    # cell_density_vector[is.na(cell_density_vector)] <- 0

    cell_abundance_val <- terra::values(N_j_raster, mat = FALSE)
    cell_abundance_val[is.na(cell_abundance_val)] <- 0
    

    obs_par_list <- as.list(current_parameter_set[, c("obs_intercept", obs_cov_names)])
    names(obs_par_list)[1] <- "intercept"

    for (sim_num in 1:n_simulations) {
      cat(sprintf("  --- Sim %d / %d ---\n", sim_num, n_simulations))

      # Retrieve the W matrix for the current reference method
      current_w_matrix <- all_w_matrices[[current_clustering_method]]
      
      # --- 6.2 Simulate Data (OPTIMIZED) ---
      train_data <- simulate_train_data(
        reference_clustering_df = current_reference_dataframe,
        obs_cov_names = obs_cov_names,
        obs_par_list = obs_par_list,
        w_matrix = current_w_matrix,           # Pass Matrix
        cell_density_vector = cell_abundance_val # Pass Density
      )
      
      # Test data simulation 
      test_data_full <- simulate_test_data(
        base_test_df = base_test_df,
        obs_cov_names = obs_cov_names,
        obs_par_list = obs_par_list,
        w_matrix = w_matrix_test,
        cell_density_vector = cell_abundance_val
      )

      # Dataset Stats
      ds_stats <- summarize_datasets(train_data, test_data_full)
      ds_stats$reference_method <- current_clustering_method
      ds_stats$param_set <- param_idx
      ds_stats$sim_num <- sim_num
      all_dataset_stats[[length(all_dataset_stats) + 1]] <- ds_stats
      
      # --- 6.3 Pre-generate Test Splits ---
      # Consistent testing across all methods
      test_splits_list <- list()
      for (r in 1:n_test_repeats) {
        test_splits_list[[r]] <- spatial_subsample_dataset(test_data_full, res_m/1000, r)
      }

      # --- 6.4 Method Loop (Fit ONCE) ---
      methods_to_test <- unique(c(current_clustering_method, comparison_method_list))
      
      # === 4. LOOP OVER METHODS (Fit Model ONCE per Sim) ===
      # for (method_name in methods_to_test) {
      
      #   cat(sprintf("\n  [Sim %d, Param %d] Fitting Method: %s (Ref: %s)\n", 
      #         sim_num, param_idx, method_name, current_clustering_method))

      #   # === 4.1. PREPARE occuN DATA ===
      #   current_clustering_df <- all_clusterings[[method_name]]
      #   if (is.list(current_clustering_df) && "result_df" %in% names(current_clustering_df)) {
      #      current_clustering_df <- current_clustering_df$result_df
      #   }
        
      #   w_matrix <- all_w_matrices[[method_name]]
        
      #   if (is.null(w_matrix)) {
      #     cat(sprintf("    Skipping %s (No W matrix)\n", method_name)); next
      #   }

      #   # Use modular function
      #   umf <- prepare_occuN_data(train_data, current_clustering_df, w_matrix, obs_cov_names, full_raster_covs)

      #   # === 4.2. FIT MODEL ===
      #   cat(sprintf("    Fitting %s (M=%d)... ", method_name, nrow(umf@y)))
        
      #   obs_formula <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
      #   state_formula <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))
        
      #   # Use modular function
      #   fm <- fit_occuN_model(umf, state_formula, obs_formula, n_reps = n_fit_repeats, stable_reps = n_fit_repeats, optimizer = selected_optimizer)
        
      #   # Define desired column names for coefficients
      #   state_col_names <- c("state_intercept", state_cov_names)
      #   obs_col_names   <- c("obs_intercept", obs_cov_names)
        
      #   # === 4.3. HANDLE FAILURE ===
      #   if (is.null(fm)) {
      #     cat("FAILED.\n")
          
      #     # Record FAILURE in Parameters
      #     na_param_row <- data.frame(
      #       reference_method = current_clustering_method, 
      #       param_set = param_idx, 
      #       sim_num = sim_num,
      #       comparison_method = method_name, 
      #       nll = NA, 
      #       convergence = 1 # Mark as failed
      #     )
      #     # Fill coefs with NA
      #     for(col in c(state_col_names, obs_col_names)) na_param_row[[col]] <- NA
      #     all_param_results[[length(all_param_results) + 1]] <- na_param_row
          
      #     # Record FAILURE in Predictions (Optional: or just skip)
      #     for(r in 1:n_test_repeats) {
      #       na_pred_row <- data.frame(
      #         reference_method = current_clustering_method,
      #         param_set = param_idx,
      #         sim_num = sim_num,
      #         comparison_method = method_name,
      #         test_repeat = r,
      #         auc = NA,
      #         auprc = NA
      #       )
      #       all_pred_results[[length(all_pred_results) + 1]] <- na_pred_row
      #     }
      #     next
      #   }
      #   cat("Done.\n")

      #   # === 4.4. SAVE PARAMETERS (Run ONCE per fit) ===
      #   est_alphas <- coef(fm, 'det')   
      #   est_betas <- coef(fm, 'state')  

      #   param_row <- data.frame(
      #     reference_method = current_clustering_method,
      #     param_set = param_idx,
      #     sim_num = sim_num,
      #     comparison_method = method_name,
      #     nll = fm@negLogLike,
      #     convergence = 0
      #   )
        
      #   # Store State Coefficients
      #   for(i in seq_along(state_col_names)) param_row[[state_col_names[i]]] <- est_betas[i]
      #   # Store Obs Coefficients
      #   for(i in seq_along(obs_col_names)) param_row[[obs_col_names[i]]] <- est_alphas[i]
        
      #   all_param_results[[length(all_param_results) + 1]] <- param_row

        
      #   # === 4.5. PREDICT & TEST (Repeat N times) ===
      #   for (repeat_num in 1:n_test_repeats) {
      #     test_df <- test_splits_list[[repeat_num]]
          
      #     # Prediction Logic
      #     X_state <- model.matrix(state_formula, data = test_df)
      #     X_obs <- model.matrix(obs_formula, data = test_df)
          
      #     pred_psi <- 1 - exp(-(exp(X_state %*% est_betas) * test_df$area_j))
      #     pred_det <- plogis(X_obs %*% est_alphas)
      #     pred_obs_prob <- pred_psi * pred_det 
          
      #     # Metrics
      #     metrics <- calculate_classification_metrics(pred_obs_prob, test_df$species_observed)
          
      #     # Store Prediction Info Only
      #     pred_row <- data.frame(
      #       reference_method = current_clustering_method,
      #       param_set = param_idx,
      #       sim_num = sim_num,
      #       comparison_method = method_name,
      #       test_repeat = repeat_num, # Key differentiator
      #       auc = metrics$auc,
      #       auprc = metrics$auprc
      #     )
          
      #     all_pred_results[[length(all_pred_results) + 1]] <- pred_row
      #   }

      #   # 6. CRITICAL MEMORY CLEANUP
      #   rm(umf, fm) # Remove heavy objects
      #   gc() # Force garbage collection
        
      # } # End Method Loop

      # Periodic cleanup after every simulation
      rm(train_data, test_data_full, test_splits_list)
      gc()
    } # End Sim Loop
  } # End Param Loop
} # End Clustering Loop

# Save Final Results
write.csv(dplyr::bind_rows(all_dataset_stats), file.path(output_dir, "dataset_descriptive_stats.csv"), row.names = FALSE)

# Save Parameters
write.csv(dplyr::bind_rows(all_param_results), file.path(output_dir, "estimated_parameters.csv"), row.names = FALSE)

# Save Predictions
write.csv(dplyr::bind_rows(all_pred_results), file.path(output_dir, "predictive_performance.csv"), row.names = FALSE)

cat("Done.")