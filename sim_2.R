# -----------------------------------------------------------------
# sim_2.R
# Complexity Gradient: Bridging Simulation and Reality
# 4 Variants with full Stats, Plotting, and 2017/2018 Real Data Split
# Uniform Sampling matches Real Data N and respects Boundary Shapefile
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
library(terra)
library(sf)
library(Matrix)
library(PRROC)

# Source existing helpers
source(file.path("R", "utils.R"))
source(file.path("R", "simulation_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))
source(file.path("R", "plotting_helpers.R"))

set.seed(123)

###
# 2. CONFIGS & VARIANTS
###

# Load Parameters and Clusterings
sim_params <- read.csv(file.path("config", "simulation_parameters.csv"))
sim_clusterings <- read.csv(file.path("config", "simulation_clusterings.csv"))

# Methods to run
comparison_method_list <- c() 
reference_method_list <- unique(sim_clusterings$method)
all_method_names <- unique(c(reference_method_list, comparison_method_list))

# Define the 4 Variants
variants <- list(
  "V1_Uniform_Simple" = list(
    loc_type = "uniform",
    state_covs_used = c("elevation"), 
    n_obs_covs = 1
  ),
  "V2_Uniform_ComplexState" = list(
    loc_type = "uniform",
    state_covs_used = c("elevation", "TCB", "TCG", "TCW", "TCA"),
    n_obs_covs = 1
  ),
  "V3_Uniform_ComplexAll" = list(
    loc_type = "uniform",
    state_covs_used = c("elevation", "TCB", "TCG", "TCW", "TCA"),
    n_obs_covs = 5
  ),
  "V4_RealLocs_ComplexAll" = list(
    loc_type = "real",
    state_covs_used = c("elevation", "TCB", "TCG", "TCW", "TCA"),
    n_obs_covs = 5
  )
)

# Global Simulation Settings
n_simulations <- 10  
n_fit_repeats <- 10
n_test_repeats <- 10

n_simulations <- 1  
n_fit_repeats <- 30
n_test_repeats <- 1

selected_optimizer <- "nlminb"
buffer_m <- 200
res_m <- 100
PARAM_LOWER <- -10
PARAM_UPPER <- 10

# Main Output Directory
main_output_dir <- file.path("simulation_experiments", "output", "sim_2")
if (!dir.exists(main_output_dir)) dir.create(main_output_dir, recursive = TRUE)

###
# 3. PREPROCESS RASTER & BOUNDARY (Global)
###

# Standardize Raster Once
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# Load Native Raster
state_cov_names_all <- names(sim_params)[2:6] 
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster_raw) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster_raw) <- state_cov_names_all

# Project to Albers (Raw Values)
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = res_m)

# Standardize
standardization_results <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- standardization_results$raster
state_cov_params <- standardization_results$params

# Area Raster
cell_area_km2 <- (res_m / 1000) * (res_m / 1000)
area_j_raster <- cov_tif_albers[[1]] * 0 + cell_area_km2
names(area_j_raster) <- "area"

# Full Covariate Dataframe (for fitting)
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))
full_raster_covs[is.na(full_raster_covs)] <- 0

# Boundary for Sampling & Plotting
boundary_shapefile_path <- file.path("state_covariate_raster", "boundary", "boundary.shp")
boundary_vect <- terra::vect(boundary_shapefile_path)
boundary_vect_albers <- terra::project(boundary_vect, albers_crs_str)

###
# 4. HELPER: DATA GENERATION
###

generate_variant_base_data <- function(variant_cfg, cov_raster_obj, boundary_vect_proj, real_train_file, real_test_file, is_test=FALSE) {
  
  # 1. Determine Target N based on Real Files
  target_file <- if(is_test) real_test_file else real_train_file
  real_df <- read.csv(target_file)
  real_df <- real_df[!is.na(real_df$latitude) & !is.na(real_df$longitude),]
  
  target_n <- nrow(real_df)
  
  if (variant_cfg$loc_type == "uniform") {
    
    cat(sprintf("    [DataGen] Uniform Sampling %d locations within Boundary...\n", target_n))
    
    # Mask the raster with the boundary to ensure we only sample valid cells INSIDE the boundary
    # (Assuming cov_raster_obj is already in same CRS as boundary_vect_proj)
    masked_raster <- terra::mask(cov_raster_obj[[1]], boundary_vect_proj)
    
    # Get valid cells from masked raster
    valid_cells <- terra::cells(masked_raster)
    
    if(length(valid_cells) < target_n) {
        warning("    [DataGen] Warning: Target N exceeds valid cells in boundary. Sampling with replacement.")
    }
    
    sampled_indices <- sample(valid_cells, target_n, replace = TRUE) 
    
    # Get center coordinates
    coords_proj <- terra::xyFromCell(cov_raster_obj[[1]], sampled_indices)
    
    # Add Jitter (+/- 20% of resolution)
    r_res <- terra::res(cov_raster_obj)[1]
    jitter_amount <- r_res * 0.2 
    coords_proj[,1] <- coords_proj[,1] + runif(nrow(coords_proj), -jitter_amount, jitter_amount)
    coords_proj[,2] <- coords_proj[,2] + runif(nrow(coords_proj), -jitter_amount, jitter_amount)
    
    df_proj <- data.frame(x = coords_proj[,1], y = coords_proj[,2])
    v_proj <- terra::vect(df_proj, geom=c("x", "y"), crs = terra::crs(cov_raster_obj))
    v_geo <- terra::project(v_proj, "+proj=longlat +datum=WGS84")
    coords_geo <- terra::crds(v_geo)
    
    prefix <- if(is_test) "test_unif" else "train_unif"
    
    base_df <- data.frame(
      checklist_id = paste0(prefix, "_", 1:target_n),
      locality_id = paste0("loc_", prefix, "_", 1:target_n),
      latitude = coords_geo[,2],
      longitude = coords_geo[,1],
      observation_date = if(is_test) "2018-06-01" else "2017-06-01",
      formatted_date = if(is_test) "2018-06-01" else "2017-06-01"
    )
    
    # Extract State Covs
    env_df <- extract_state_covs(base_df, cov_raster_obj)
    base_df <- dplyr::inner_join(base_df, env_df, by = "checklist_id")
    
  } else {
    # REAL LOCATIONS
    cat(sprintf("    [DataGen] Using Real Data (N=%d)...\n", target_n))
    
    base_df <- real_df[, c("checklist_id", "locality_id", "latitude", "longitude", "observation_date")]
    base_df$formatted_date <- base_df$observation_date
    
    # Extract State Covs
    env_df <- extract_state_covs(base_df, cov_raster_obj)
    base_df <- dplyr::inner_join(base_df, env_df, by = "checklist_id")
  }
  
  # 2. Simulate Observation Covariates (Placeholders)
  obs_cov_names <- c("duration_minutes", "effort_distance_km", "number_observers", "time_observations_started", "day_of_year")
  for(col in obs_cov_names){
    base_df[[col]] <- rnorm(nrow(base_df)) 
  }
  
  return(base_df)
}

###
# 5. MAIN VARIANT LOOP
###

# File paths for Real Data Counts
file_train_2017 <- file.path("checklist_data", "species", "AMCR", "AMCR_zf_filtered_region_2017.csv")
file_test_2018  <- file.path("checklist_data", "species", "AMCR", "AMCR_zf_filtered_region_2018.csv")

all_param_results <- list()
all_pred_results <- list()

# Define plot order
all_method_names_plot_order <- c(
  "0.25-kmSq", "0.5-kmSq", "1-kmSq", "2-kmSq",
  "clustGeo-50-80", "clustGeo-50-40", "clustGeo-50-20", "clustGeo-50-10", "DBSC"
)

for (v_name in names(variants)) {
  
  variant <- variants[[v_name]]
  
  # Variant Output Directory
  v_output_dir <- file.path(main_output_dir, v_name)
  if (!dir.exists(v_output_dir)) dir.create(v_output_dir, recursive = TRUE)
  
  cat(paste("\n################################################\n"))
  cat(paste("### STARTING VARIANT:", v_name, "###\n"))
  cat(paste("################################################\n"))
  
  current_state_covs <- variant$state_covs_used
  all_obs_headers <- c("duration_minutes", "effort_distance_km", "number_observers", "time_observations_started", "day_of_year")
  current_obs_covs <- all_obs_headers[1:variant$n_obs_covs] 
  
  # --- A. Generate Base Data ---
  # Passing file paths so generator can match N exactly
  cat("  --- Generating Variant Dataset ---\n")
  base_train_df <- generate_variant_base_data(
      variant, cov_tif_albers, boundary_vect_albers, 
      file_train_2017, file_test_2018, is_test=FALSE
  )
  base_test_df  <- generate_variant_base_data(
      variant, cov_tif_albers, boundary_vect_albers, 
      file_train_2017, file_test_2018, is_test=TRUE
  )
  
  # --- B. Pre-compute Clusterings (Train) ---
  cat("  --- Pre-computing Clusterings (Train) ---\n")
  all_clusterings <- get_clusterings(all_method_names, base_train_df, names(state_cov_raster_raw), NULL)
  
  # --- C. Initial Geometries (No Splitting Yet) ---
  cat("  --- Creating Initial Geometries ---\n")
  all_site_geometries <- list()
  for (method_name in all_method_names) {
    cluster_data <- all_clusterings[[method_name]]
    if (is.list(cluster_data) && "result_df" %in% names(cluster_data)) cluster_data <- cluster_data$result_df
    
    if(!is.null(cluster_data)){
       all_site_geometries[[method_name]] <- create_site_geometries(cluster_data, cov_tif_albers, buffer_m, method_name)
    }
  }
  
  # --- D. Stats & Plots: PRE-SPLIT ---
  cat("  --- Stats & Plots: PRE-SPLIT ---\n")
  
  # 1. Stats
  summary_pre <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
  write.csv(summary_pre, file.path(v_output_dir, "clustering_stats_PRE_split.csv"), row.names = FALSE)
  
  similarity_pre <- compute_clustering_similarity(all_clusterings, reference_method_list, all_method_names)
  write.csv(similarity_pre, file.path(v_output_dir, "similarity_stats_PRE_split.csv"), row.names = FALSE)
  
  # 2. Plots
  plot_sites(
    base_train_df = base_train_df,
    all_clusterings = all_clusterings,
    all_site_geometries = all_site_geometries,
    elevation_raster = cov_tif_albers_raw, 
    methods_to_plot = all_method_names_plot_order,
    boundary_shp_path = boundary_shapefile_path,
    output_path = file.path(v_output_dir, "visualization_PRE.png"),
    cluster_labels = TRUE
  )

  # --- E. Split Disjoint Sites & W Matrices ---
  cat("  --- Splitting Disjoint Sites & Generating W ---\n")
  all_w_matrices <- list()
  
  for (method_name in names(all_site_geometries)) {
    curr_geoms <- all_site_geometries[[method_name]]
    curr_data_obj <- all_clusterings[[method_name]]
    is_list <- is.list(curr_data_obj) && "result_df" %in% names(curr_data_obj)
    curr_data <- if(is_list) curr_data_obj$result_df else curr_data_obj
    
    # Split
    split_res <- disjoint_site_geometries(curr_geoms, curr_data)
    
    # Update
    all_site_geometries[[method_name]] <- split_res$geoms
    if (is_list) {
      all_clusterings[[method_name]]$result_df <- split_res$data
    } else {
      all_clusterings[[method_name]] <- split_res$data
    }
    
    # W Matrix
    all_w_matrices[[method_name]] <- generate_overlap_matrix(split_res$geoms, cov_tif_albers)
  }

  # --- F. Stats & Plots: POST-SPLIT ---
  cat("  --- Stats & Plots: POST-SPLIT ---\n")
  
  # 1. Stats
  summary_post <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
  write.csv(summary_post, file.path(v_output_dir, "clustering_stats_POST_split.csv"), row.names = FALSE)
  
  similarity_post <- compute_clustering_similarity(all_clusterings, reference_method_list, all_method_names)
  write.csv(similarity_post, file.path(v_output_dir, "similarity_stats_POST_split.csv"), row.names = FALSE)
  
  # 2. Plots
  plot_sites(
    base_train_df = base_train_df,
    all_clusterings = all_clusterings,
    all_site_geometries = all_site_geometries,
    elevation_raster = cov_tif_albers_raw, 
    methods_to_plot = all_method_names_plot_order,
    boundary_shp_path = boundary_shapefile_path,
    output_path = file.path(v_output_dir, "visualization_POST.png"),
    cluster_labels = TRUE
  )
  
  # --- G. Clean up Raster Memory ---
  rm(all_site_geometries)
  gc()

  # --- H. Prepare Test Structures ---
  test_structures <- prepare_test_spatial_structures(
    test_df = base_test_df,
    albers_crs = albers_crs_str,
    buffer_m = buffer_m,
    cov_raster_albers = cov_tif_albers,
    area_raster = area_j_raster
  )
  base_test_df_ready <- test_structures$test_df
  w_matrix_test <- test_structures$w_matrix
  rm(test_structures)
  
  
  # --- I. Simulation Loops ---
  for (ref_method in reference_method_list) {
    cat(paste("\n  === REF CLUSTERING:", ref_method, "===\n"))
    
    current_reference_dataframe <- all_clusterings[[ref_method]]
    if (is.list(current_reference_dataframe)) current_reference_dataframe <- current_reference_dataframe$result_df
    
    current_w_matrix <- all_w_matrices[[ref_method]]

    for (param_idx in seq_len(nrow(sim_params))) {
      curr_params <- sim_params[param_idx, ]
      
      # 1. Extract Variant-Specific State Params
      state_par_list <- list(intercept = curr_params$state_intercept)
      for(sc in names(full_raster_covs)){
        if(sc %in% current_state_covs){
          state_par_list[[sc]] <- curr_params[[sc]]
        } else {
          state_par_list[[sc]] <- 0 
        }
      }
      
      # 2. Density Surface
      log_lambda_j <- cov_tif_albers[[1]] * 0 + state_par_list$intercept
      for (nm in names(full_raster_covs)) {
        if(state_par_list[[nm]] != 0) {
           log_lambda_j <- log_lambda_j + (cov_tif_albers[[nm]] * state_par_list[[nm]])
        }
      }
      cell_density_val <- terra::values(exp(log_lambda_j), mat = FALSE)
      cell_density_val[is.na(cell_density_val)] <- 0
      
      # 3. Extract Variant-Specific Obs Params
      obs_par_list <- list(intercept = curr_params$obs_intercept)
      for(oc in current_obs_covs){
        obs_par_list[[oc]] <- curr_params[[oc]]
      }

      for (sim_num in 1:n_simulations) {
        set.seed(sim_num * param_idx) 
        
        # Simulate Data
        train_data <- simulate_train_data(
          reference_clustering_df = current_reference_dataframe,
          obs_cov_names = current_obs_covs, 
          obs_par_list = obs_par_list,
          w_matrix = current_w_matrix,
          cell_density_vector = cell_density_val
        )
        
        test_data_full <- simulate_test_data(
          base_test_df = base_test_df_ready,
          obs_cov_names = current_obs_covs,
          obs_par_list = obs_par_list,
          w_matrix = w_matrix_test,
          cell_density_vector = cell_density_val
        )
        
        # Test Splits
        test_splits_list <- list()
        for (r in 1:n_test_repeats) {
          test_splits_list[[r]] <- spatial_subsample_dataset(test_data_full, res_m/1000, r)
        }
        
        methods_to_test <- unique(c(ref_method, comparison_method_list))
        
        for (method_name in methods_to_test) {
          
          fit_clustering_df <- all_clusterings[[method_name]]
          if (is.list(fit_clustering_df)) fit_clustering_df <- fit_clustering_df$result_df
          fit_w_matrix <- all_w_matrices[[method_name]]
          
          if(is.null(fit_w_matrix)) next
          
          umf <- prepare_occuN_data(train_data, fit_clustering_df, fit_w_matrix, current_obs_covs, full_raster_covs)
          
          state_form <- as.formula(paste("~", paste(current_state_covs, collapse = " + ")))
          obs_form <- as.formula(paste("~", paste(current_obs_covs, collapse = " + ")))
          
          fm <- fit_occuN_model(
            umf, state_form, obs_form, 
            n_reps = n_fit_repeats, stable_reps = 3, 
            optimizer = selected_optimizer,
            lower = PARAM_LOWER, upper = PARAM_UPPER
          )
          
          if (!is.null(fm)) {
            est_betas <- coef(fm, 'state')
            est_alphas <- coef(fm, 'det')
            
            # Save Params
            param_row <- data.frame(
              Variant = v_name,
              ref_method = ref_method,
              comp_method = method_name,
              param_set = param_idx,
              sim_num = sim_num,
              NLL = fm@negLogLike,
              convergence = 0
            )
            for(n in names(est_betas)) param_row[[paste0("state_", n)]] <- est_betas[[n]]
            for(n in names(est_alphas)) param_row[[paste0("obs_", n)]] <- est_alphas[[n]]
            
            all_param_results[[length(all_param_results) + 1]] <- param_row
            
            # Predict & Test
            for (r in 1:n_test_repeats) {
              test_df <- test_splits_list[[r]]
              
              X_state <- model.matrix(state_form, data = test_df)
              X_obs <- model.matrix(obs_form, data = test_df)
              
              if(ncol(X_state) != length(est_betas)) next 
              
              pred_psi <- 1 - exp(-(exp(X_state %*% est_betas) * test_df$area_j))
              pred_det <- plogis(X_obs %*% est_alphas)
              pred_obs_prob <- pred_psi * pred_det 
              
              metrics <- calculate_classification_metrics(pred_obs_prob, test_df$species_observed)
              
              pred_row <- data.frame(
                Variant = v_name,
                ref_method = ref_method,
                comp_method = method_name,
                param_set = param_idx,
                sim_num = sim_num,
                test_repeat = r,
                auc = metrics$auc,
                auprc = metrics$auprc
              )
              all_pred_results[[length(all_pred_results) + 1]] <- pred_row
            }
          }
          rm(fm, umf); gc()
        } 
        rm(train_data, test_data_full, test_splits_list); gc()
      } 
    } 
  } 
  
  rm(all_clusterings, all_w_matrices, base_train_df, base_test_df_ready); gc()
  
} # End Variant Loop

###
# 6. SAVE GLOBAL RESULTS
###

final_params <- dplyr::bind_rows(all_param_results)
final_preds <- dplyr::bind_rows(all_pred_results)

write.csv(final_params, file.path(main_output_dir, "gradient_parameters.csv"), row.names = FALSE)
write.csv(final_preds, file.path(main_output_dir, "gradient_predictions.csv"), row.names = FALSE)

cat("\n--- Gradient Experiment Complete ---\n")