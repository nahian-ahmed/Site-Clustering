# -----------------------------------------------------------------
# Simulation for occuN model (Optimized & Parallelized)
# -----------------------------------------------------------------

# --- 1. Setup & Dependencies ---
install_now = FALSE # Set to TRUE if first run
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
library(Matrix)
library(doParallel) # For parallel processing
library(foreach)    # For parallel loops
library(dggridR)    # Ensure this is loaded

# Load Helpers
source(file.path("R", "utils.R"))
source(file.path("R", "simulation_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))
source(file.path("R", "plotting_helpers.R"))

set.seed(123) 

# --- 2. Comparison Methods Configuration ---
comparison_method_list <- c(
  "1to10", "2to10", "2to10-sameObs", "1-kmSq",
  "lat-long", "rounded-4", "SVS", "1-per-UL",
  "DBSC", "BayesOptClustGeo"
)

selected_optimizer <- "nlminb"

# --- 3. Load Simulation Config ---
sim_params <- read.delim(file.path("config", "simulation_parameters.csv"), sep = ",", header = T)
sim_clusterings <- read.delim(file.path("config", "simulation_clusterings.csv"), sep = ",", header = T)

# Simulation Settings
n_simulations <- 25
n_fit_repeats <- 10 # Reduced from 25 (Early stopping handles the rest)
n_test_repeats <- 25

res_m <- 30 
buffer_m <- 150 

state_cov_names <- names(sim_params)[2:6]
obs_cov_names <- names(sim_params)[8:12]

# --- 4. Pre-processing (Static Data) ---
state_cov_raster <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster) <- state_cov_names

boundary_shapefile_path <- file.path("state_covariate_raster", "boundary", "boundary.shp")

base_train_data <- prepare_train_data(state_cov_names, obs_cov_names, state_cov_raster)
base_train_df <- base_train_data$train_df
norm_list <- base_train_data$norm_list

base_test_df <- prepare_test_data(state_cov_names, obs_cov_names, state_cov_raster, norm_list)

# Pre-calculate Albers info
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
cov_tif_albers <- terra::project(state_cov_raster, albers_crs_str, method="bilinear", res = res_m)
area_j_raster <- terra::cellSize(cov_tif_albers, unit="m")
full_raster_covs <- as.data.frame(terra::values(state_cov_raster))[, state_cov_names, drop = FALSE]
full_raster_covs[is.na(full_raster_covs)] <- 0

# --- 5. Pre-compute Clusterings, Geometries & WEIGHT MATRICES ---
reference_method_list <- sim_clusterings$method
all_method_names <- unique(c(reference_method_list, comparison_method_list))

cat("--- Pre-computing ALL clusterings... ---\n")
all_clusterings <- get_clusterings(all_method_names, base_train_df, state_cov_names, NULL)

cat("--- Pre-computing site geometries... ---\n")
all_site_geometries <- list() 
for (method_name in all_method_names) {
  cluster_data <- all_clusterings[[method_name]]
  if (is.list(cluster_data) && "result_df" %in% names(cluster_data)) cluster_data <- cluster_data$result_df
  
  if (!is.null(cluster_data)) {
    all_site_geometries[[method_name]] <- create_site_geometries(cluster_data, state_cov_raster, buffer_m, method_name, "km")
  }
}

# --- OPTIMIZATION: Pre-compute Sparse Weight Matrices ---
# This replaces repeated spatial extraction in the loop
cat("--- Pre-computing Spatial Weight Matrices (Optimization) ---\n")
site_weight_matrices <- list()
dummy_raster <- cov_tif_albers[[1]] # Used for cell indexing

for (method_name in names(all_site_geometries)) {
  geoms <- all_site_geometries[[method_name]]
  
  # Extract cell numbers and overlap fractions ONCE
  extraction <- terra::extract(dummy_raster, geoms, cells=TRUE, exact=TRUE, ID=TRUE)
  
  if(nrow(extraction) > 0) {
    # Create sparse matrix: Rows = Sites, Cols = Raster Cells
    site_weight_matrices[[method_name]] <- Matrix::sparseMatrix(
      i = extraction$ID,
      j = extraction$cell,
      x = extraction$fraction, 
      dims = c(nrow(geoms), terra::ncell(dummy_raster)),
      # +++ CRITICAL FIX: Add Row Names (Site IDs) +++
      dimnames = list(as.character(geoms$site), NULL)
    )
  }
}

# Save clustering stats
clustering_summary_df <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.csv(clustering_summary_df, file.path(output_dir, "clustering_descriptive_stats.csv"), row.names = FALSE)

# === CALCULATE CLUSTERING SIMILARITY STATISTICS ===
cat("--- Calculating pairwise clustering similarity metrics (ARI, AMI, NID) ---\n")

clustering_similarity_list <- list()

# We want to compare every Reference method against every Comparison method
# (and potentially Reference vs Reference if you want to check stability)

# 1. Define the pairs to compare
# Usually: Reference methods (Sim Truth) vs All other methods
ref_methods_to_check <- sim_clusterings$method
comp_methods_to_check <- unique(c(sim_clusterings$method, comparison_method_list))

for (ref_name in ref_methods_to_check) {
  for (comp_name in comp_methods_to_check) {
  
  # Extract Dataframes
  ref_dat <- all_clusterings[[ref_name]]
  if (is.list(ref_dat) && "result_df" %in% names(ref_dat)) ref_dat <- ref_dat$result_df
  
  comp_dat <- all_clusterings[[comp_name]]
  if (is.list(comp_dat) && "result_df" %in% names(comp_dat)) comp_dat <- comp_dat$result_df
  
  if (is.null(ref_dat) || is.null(comp_dat)) next
  
  # Calculate Stats
  stats <- calculate_clustering_stats(ref_dat, comp_dat)
  
  # Store
  clustering_similarity_list[[length(clustering_similarity_list) + 1]] <- data.frame(
    reference_method = ref_name,
    comparison_method = comp_name,
    ARI = stats$ARI,
    AMI = stats$AMI,
    NID = stats$NID
  )
  }
}

# Bind and Save
clustering_similarity_df <- dplyr::bind_rows(clustering_similarity_list)
write.csv(clustering_similarity_df, file.path(output_dir, "clustering_similarity_stats.csv"), row.names = FALSE)
cat(sprintf("--- Clustering similarity stats saved to %s/clustering_similarity_stats.csv ---\n", output_dir))


all_method_names_plot_order <- c(
  "1to10", "2to10", "2to10-sameObs", "lat-long", "SVS", "1-per-UL", 
  "2-kmSq", "1-kmSq", "0.5-kmSq", "0.25-kmSq", "0.125-kmSq", "rounded-4",
  "clustGeo-50-20", "clustGeo-50-40", "clustGeo-50-60", "clustGeo-50-80", "BayesOptClustGeo", "DBSC"
)
# Call plot function
site_plot <- plot_sites(
  base_train_df = base_train_df,
  all_clusterings = all_clusterings,
  all_site_geometries = all_site_geometries,
  elevation_raster = state_cov_raster,
  methods_to_plot = all_method_names_plot_order,
  boundary_shp_path = boundary_shapefile_path,
  output_path = file.path(output_dir, "site_cluster_visualization.png")
)


# --- 6. Main Simulation Loop (PARALLELIZED) ---

# Setup Parallel Backend
num_cores <- parallel::detectCores() - 1
registerDoParallel(cores = num_cores)
cat(sprintf("--- Starting Simulation on %d cores ---\n", num_cores))

# Storage lists (will accumulate results after loops)
all_dataset_stats <- list()
all_param_results <- list() 
all_pred_results <- list()

for (cluster_idx in seq_len(nrow(sim_clusterings))) {
  current_clustering_method <- sim_clusterings$method[cluster_idx]
  cat(paste("\n=== REF CLUSTERING:", current_clustering_method, "===\n"))
  
  current_reference_dataframe <- all_clusterings[[current_clustering_method]]
  current_site_geometries <- all_site_geometries[[current_clustering_method]]
  
  # Get the pre-calculated weight matrix for the REFERENCE method
  current_ref_weight_matrix <- site_weight_matrices[[current_clustering_method]]

  for (param_idx in seq_len(nrow(sim_params))) {
    current_parameter_set <- sim_params[param_idx, ]
    cat(paste("  Param Set:", param_idx, "\n"))

    # --- 6.1 Pre-calc N_j_raster ---
    state_par_list <- as.list(current_parameter_set[, c("state_intercept", state_cov_names)])
    names(state_par_list)[1] <- "intercept"
    
    log_lambda_j <- cov_tif_albers[[1]] * 0 + state_par_list$intercept
    for (nm in state_cov_names) log_lambda_j <- log_lambda_j + (cov_tif_albers[[nm]] * state_par_list[[nm]])
    
    N_j_raster <- exp(log_lambda_j) * area_j_raster
    
    obs_par_list <- as.list(current_parameter_set[, c("obs_intercept", obs_cov_names)])
    names(obs_par_list)[1] <- "intercept"
    
    # Wrap rasters for passing to parallel workers
    packed_N_j <- terra::wrap(N_j_raster)
    packed_area_j <- terra::wrap(area_j_raster)

    # --- 6.2 Parallel Simulation Loop ---
    # We use foreach to distribute the 'sim_num' loop
    sim_results_list <- foreach(
      sim_num = 1:n_simulations,
      .packages = c('terra', 'dplyr', 'unmarked', 'Matrix', 'rje', 'PRROC', 'dggridR'),
      .export = c('simulate_train_data', 'simulate_test_data', 
                  'prepare_occuN_data', 'fit_occuN_model',
                  'calculate_weighted_sum', 'norm_ds', 
                  'spatial_subsample_dataset', 'summarize_datasets',
                  'calculate_classification_metrics') 
    ) %dopar% {
      
      # Unwrap rasters inside the worker
      N_j_worker <- terra::unwrap(packed_N_j)
      area_j_worker <- terra::unwrap(packed_area_j)
      
      # Local storage for this specific simulation run
      local_stats <- list()
      local_params <- list()
      local_preds <- list()
      
      # --- A. Simulate Data ---
      # Use the Optimized simulate_train_data with weight_matrix
      train_data <- simulate_train_data(
        reference_clustering_df = current_reference_dataframe,
        site_geoms_sf = current_site_geometries,
        obs_cov_names = obs_cov_names,
        obs_par_list = obs_par_list,
        N_j_raster = N_j_worker,
        weight_matrix = current_ref_weight_matrix # <--- OPTIMIZATION
      )
      
      test_data_full <- simulate_test_data(
        base_test_df = base_test_df,
        obs_cov_names = obs_cov_names,
        obs_par_list = obs_par_list,
        N_j_raster = N_j_worker,
        albers_crs_str = albers_crs_str,
        area_j_raster = area_j_worker
      )
      
      # Stats
      ds_stats <- summarize_datasets(train_data, test_data_full)
      ds_stats$reference_method <- current_clustering_method
      ds_stats$param_set <- param_idx
      ds_stats$sim_num <- sim_num
      local_stats[[1]] <- ds_stats
      
      # --- B. Test Splits ---
      test_splits_list <- list()
      for (r in 1:n_test_repeats) {
        test_splits_list[[r]] <- spatial_subsample_dataset(test_data_full, res_m/1000, r)
      }
      
      # --- C. Method Loop (Fit & Predict) ---
      methods_to_test <- unique(c(current_clustering_method, comparison_method_list))
      
      for (method_name in methods_to_test) {
        
        # 1. Prepare Data
        current_clustering_df <- all_clusterings[[method_name]]
        if (is.list(current_clustering_df) && "result_df" %in% names(current_clustering_df)) {
          current_clustering_df <- current_clustering_df$result_df
        }
        
        # Get geometries and Weight Matrix for the *fitting* method
        current_geoms <- all_site_geometries[[method_name]]
        # W matrix is now stored in our pre-calculated list, NOT as an attribute
        w_matrix <- site_weight_matrices[[method_name]] 
        
        if (is.null(w_matrix)) next
        
        umf <- prepare_occuN_data(train_data, current_clustering_df, w_matrix, obs_cov_names, full_raster_covs)
        
        obs_formula <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
        state_formula <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))
        
        # 2. Fit Model (uses optimized fit_occuN_model with early stopping)
        fm <- fit_occuN_model(umf, state_formula, obs_formula, n_reps = n_fit_repeats, optimizer = selected_optimizer)
        
        state_col_names <- c("state_intercept", state_cov_names)
        obs_col_names   <- c("obs_intercept", obs_cov_names)
        
        # 3. Save Parameters
        param_row <- data.frame(
          reference_method = current_clustering_method,
          param_set = param_idx,
          sim_num = sim_num,
          comparison_method = method_name,
          nll = if(is.null(fm)) NA else fm@negLogLike,
          convergence = if(is.null(fm)) 1 else 0
        )
        
        if (!is.null(fm)) {
          est_alphas <- coef(fm, 'det')
          est_betas <- coef(fm, 'state')
          for(i in seq_along(state_col_names)) param_row[[state_col_names[i]]] <- est_betas[i]
          for(i in seq_along(obs_col_names)) param_row[[obs_col_names[i]]] <- est_alphas[i]
        } else {
          for(col in c(state_col_names, obs_col_names)) param_row[[col]] <- NA
        }
        local_params[[length(local_params) + 1]] <- param_row
        
        # 4. Predict & Test (only if fit succeeded)
        if (!is.null(fm)) {
          for (repeat_num in 1:n_test_repeats) {
            test_df <- test_splits_list[[repeat_num]]
            
            X_state <- model.matrix(state_formula, data = test_df)
            X_obs <- model.matrix(obs_formula, data = test_df)
            
            pred_psi <- 1 - exp(-(exp(X_state %*% est_betas) * test_df$area_j))
            pred_det <- plogis(X_obs %*% est_alphas)
            pred_obs_prob <- pred_psi * pred_det 
            
            metrics <- calculate_classification_metrics(pred_obs_prob, test_df$species_observed)
            
            pred_row <- data.frame(
              reference_method = current_clustering_method,
              param_set = param_idx,
              sim_num = sim_num,
              comparison_method = method_name,
              test_repeat = repeat_num,
              auc = metrics$auc,
              auprc = metrics$auprc
            )
            local_preds[[length(local_preds) + 1]] <- pred_row
          }
        }
      } # End Method Loop
      
      # Return lists for this sim
      return(list(
        stats = dplyr::bind_rows(local_stats),
        params = dplyr::bind_rows(local_params),
        preds = dplyr::bind_rows(local_preds)
      ))
      
    } # End Parallel Loop
    
    # --- 6.3 Aggregate Parallel Results ---
    cat("  Aggregating results...\n")
    for (res in sim_results_list) {
      if (!is.null(res$stats)) all_dataset_stats[[length(all_dataset_stats) + 1]] <- res$stats
      if (!is.null(res$params)) all_param_results[[length(all_param_results) + 1]] <- res$params
      if (!is.null(res$preds)) all_pred_results[[length(all_pred_results) + 1]] <- res$preds
    }
    
  } # End Param Loop
} # End Clustering Loop

# Stop Parallel Cluster
stopImplicitCluster()

# --- 7. Save Results ---
output_dir <- file.path("simulation_experiments", "output")
write.csv(dplyr::bind_rows(all_dataset_stats), file.path(output_dir, "dataset_descriptive_stats.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_param_results), file.path(output_dir, "estimated_parameters.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_pred_results), file.path(output_dir, "predictive_performance.csv"), row.names = FALSE)

cat("Done.")