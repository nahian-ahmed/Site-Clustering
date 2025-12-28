# -----------------------------------------------------------------
# sim_2.R
# Complexity Gradient: Bridging Simulation and Reality
# 4 Variants with full Stats, Plotting, and 2017/2018 Real Data Split
# Consistent Uniform Locations for V1-V3 (Generated Once)
# Clean Parameter Names & Dataset Descriptive Stats
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
  )
  # ,
  # "V2_Uniform_ComplexState" = list(
  #   loc_type = "uniform",
  #   state_covs_used = c("elevation", "TCB", "TCG", "TCW", "TCA"),
  #   n_obs_covs = 1
  # ),
  # "V3_Uniform_ComplexAll" = list(
  #   loc_type = "uniform",
  #   state_covs_used = c("elevation", "TCB", "TCG", "TCW", "TCA"),
  #   n_obs_covs = 5
  # ),
  # "V4_RealLocs_ComplexAll" = list(
  #   loc_type = "real",
  #   state_covs_used = c("elevation", "TCB", "TCG", "TCW", "TCA"),
  #   n_obs_covs = 5
  # )
)

# Global Simulation Settings
n_simulations <- 3
n_fit_repeats <- 30
n_test_repeats <- 3
n_stable_repeats <- 30
selected_optimizer <- "nlminb"
buffer_m <- 200
res_m <- 100

DATA_SCALE_FACTOR <- 3

PARAM_LOWER <- -10
PARAM_UPPER <- 10
INIT_LOWER <- -2
INIT_UPPER <- 2

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
# 4. DATA PREPARATION (GLOBAL)
###

obs_cov_names <- c("duration_minutes", "effort_distance_km", "number_observers", "time_observations_started", "day_of_year")

# --- A. REAL DATA (Load & Filter) ---
cat("--- Loading Real Data (Main) ---\n")

# Train (2017)
train_file <- file.path("checklist_data", "species", "AMCR", "AMCR_zf_filtered_region_2017.csv")
train_df_raw <- read.csv(train_file)
train_df_real <- train_df_raw[!is.na(train_df_raw$duration_minutes), ]
train_df_real <- train_df_real[
  train_df_real$observation_date >= "2017-05-15" & 
  train_df_real$observation_date <= "2017-07-09", 
]
real_train_master <- train_df_real[, c("checklist_id", "locality_id", "latitude", "longitude", "observation_date")]
real_train_master$formatted_date <- real_train_master$observation_date

TARGET_N_TRAIN <- nrow(real_train_master)
cat(sprintf("  -> Real Training N: %d\n", TARGET_N_TRAIN))

# --- MODIFIED BLOCK ---
# Calculate scaled N for Uniform variants (V1-V3)
SCALED_N_TRAIN <- ceiling(TARGET_N_TRAIN * DATA_SCALE_FACTOR)

cat(sprintf("  -> Generating Uniform Training Data with N: %d (Factor: %.2f)\n", SCALED_N_TRAIN, DATA_SCALE_FACTOR))

# Use SCALED_N_TRAIN for training, but keep TARGET_N_TEST for testing
unif_train_master <- generate_uniform_master(SCALED_N_TRAIN, cov_tif_albers, boundary_vect_albers, "train_unif", "2017-06-01")
unif_test_master  <- generate_uniform_master(TARGET_N_TEST,  cov_tif_albers, boundary_vect_albers, "test_unif",  "2018-06-01")
# ----------------------


# Test (2018)
test_file <- file.path("checklist_data", "species", "AMCR", "AMCR_zf_filtered_region_2018.csv")
test_df_raw <- read.csv(test_file)
test_df_real <- test_df_raw[!is.na(test_df_raw$duration_minutes), ]
test_df_real <- test_df_real[
  test_df_real$observation_date >= "2018-05-15" & 
  test_df_real$observation_date <= "2018-07-09", 
]
real_test_master <- test_df_real[, c("checklist_id", "locality_id", "latitude", "longitude", "observation_date")]
real_test_master$formatted_date <- real_test_master$observation_date

TARGET_N_TEST <- nrow(real_test_master)
cat(sprintf("  -> Real Testing N: %d\n", TARGET_N_TEST))

# --- B. UNIFORM DATA (Generate Once) ---
cat("--- Generating Uniform Data (Main) ---\n")

generate_uniform_master <- function(n, raster_obj, boundary_v, prefix, date_str) {
  # Mask raster
  masked <- terra::mask(raster_obj[[1]], boundary_v)
  valid_cells <- terra::cells(masked)
  
  if(length(valid_cells) < n) warning("Target N > Valid Cells in Boundary")
  
  # Sample
  sampled_idx <- sample(valid_cells, n, replace = TRUE)
  coords_proj <- terra::xyFromCell(raster_obj[[1]], sampled_idx)
  
  # Jitter
  r_res <- terra::res(raster_obj)[1]
  jitter_amt <- r_res * 0.2
  coords_proj[,1] <- coords_proj[,1] + runif(n, -jitter_amt, jitter_amt)
  coords_proj[,2] <- coords_proj[,2] + runif(n, -jitter_amt, jitter_amt)
  
  # To Lat/Lon
  v_proj <- terra::vect(coords_proj, crs = terra::crs(raster_obj))
  v_geo <- terra::project(v_proj, "+proj=longlat +datum=WGS84")
  coords_geo <- terra::crds(v_geo)
  
  df <- data.frame(
    checklist_id = paste0(prefix, "_", 1:n),
    locality_id = paste0("loc_", prefix, "_", 1:n),
    latitude = coords_geo[,2],
    longitude = coords_geo[,1],
    observation_date = date_str,
    formatted_date = date_str
  )
  return(df)
}

unif_train_master <- generate_uniform_master(TARGET_N_TRAIN, cov_tif_albers, boundary_vect_albers, "train_unif", "2017-06-01")
unif_test_master  <- generate_uniform_master(TARGET_N_TEST,  cov_tif_albers, boundary_vect_albers, "test_unif",  "2018-06-01")


# --- C. ENRICH DATA (Add Obs Covs & Extract State Covs) ---
cat("--- Enriching Main Datasets ---\n")

enrich_dataset <- function(df, raster_obj, obs_names) {
  # 1. Obs Covariates (Random Standard Normal)
  for(col in obs_names) df[[col]] <- rnorm(nrow(df))
  
  # 2. State Covariates (Extract from Raster)
  env_df <- extract_state_covs(df, raster_obj)
  df <- dplyr::inner_join(df, env_df, by = "checklist_id")
  return(df)
}

real_train_master <- enrich_dataset(real_train_master, cov_tif_albers, obs_cov_names)
real_test_master  <- enrich_dataset(real_test_master,  cov_tif_albers, obs_cov_names)
unif_train_master <- enrich_dataset(unif_train_master, cov_tif_albers, obs_cov_names)
unif_test_master  <- enrich_dataset(unif_test_master,  cov_tif_albers, obs_cov_names)


###
# 5. MAIN VARIANT LOOP
###

all_param_results <- list()
all_pred_results <- list()
all_dataset_stats <- list() # New list for Dataset Stats

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
  current_obs_covs <- obs_cov_names[1:variant$n_obs_covs] 
  
  # Define Clean Column Names
  state_col_names <- c("state_intercept", current_state_covs)
  obs_col_names   <- c("obs_intercept", current_obs_covs)
  
  # --- A. Select Dataset ---
  if (variant$loc_type == "uniform") {
    cat("  --- Using Main UNIFORM Data (V1/V2/V3 Consistent) ---\n")
    base_train_df <- unif_train_master
    base_test_df  <- unif_test_master
  } else {
    cat("  --- Using Main REAL Data ---\n")
    base_train_df <- real_train_master
    base_test_df  <- real_test_master
  }
  
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
    
    # Handle list vs dataframe for Reference
    current_reference_dataframe <- all_clusterings[[ref_method]]
    if (is.list(current_reference_dataframe) && "result_df" %in% names(current_reference_dataframe)) {
      current_reference_dataframe <- current_reference_dataframe$result_df
    }
    
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
        
        # --- SAVE DATASET STATS ---
        ds_stats <- summarize_datasets(train_data, test_data_full)
        ds_stats$Variant <- v_name
        ds_stats$ref_method <- ref_method
        ds_stats$param_set <- param_idx
        ds_stats$sim_num <- sim_num
        all_dataset_stats[[length(all_dataset_stats) + 1]] <- ds_stats
        
        # Test Splits
        test_splits_list <- list()
        for (r in 1:n_test_repeats) {
          test_splits_list[[r]] <- spatial_subsample_dataset(test_data_full, res_m/1000, r)
        }
        
        methods_to_test <- unique(c(ref_method, comparison_method_list))
        
        for (method_name in methods_to_test) {
          
          # Handle list vs dataframe for Fitted Model
          fit_clustering_df <- all_clusterings[[method_name]]
          if (is.list(fit_clustering_df) && "result_df" %in% names(fit_clustering_df)) {
             fit_clustering_df <- fit_clustering_df$result_df
          }
          fit_w_matrix <- all_w_matrices[[method_name]]
          
          if(is.null(fit_w_matrix)) next
          
          umf <- prepare_occuN_data(train_data, fit_clustering_df, fit_w_matrix, current_obs_covs, full_raster_covs)
          
          state_form <- as.formula(paste("~", paste(current_state_covs, collapse = " + ")))
          obs_form <- as.formula(paste("~", paste(current_obs_covs, collapse = " + ")))
          

          # --- Fit Model with Bounds ---
          
          fm <- fit_occuN_model(
            umf, state_form, obs_form, 
            n_reps = n_fit_repeats, stable_reps = n_stable_repeats, 
            optimizer = selected_optimizer,
            lower = PARAM_LOWER,
            upper = PARAM_UPPER,
            init_lower = INIT_LOWER,
            init_upper = INIT_UPPER
          )
          
          if (!is.null(fm)) {
            est_betas <- coef(fm, 'state')
            est_alphas <- coef(fm, 'det')
            
            # Init param row
            param_row <- data.frame(
              Variant = v_name,
              ref_method = ref_method,
              comp_method = method_name,
              param_set = param_idx,
              sim_num = sim_num,
              NLL = fm@negLogLike,
              convergence = 0
            )
            
            # Save Parameters with Clean Names
            # 1. State
            if(length(est_betas) == length(state_col_names)){
               for(i in seq_along(state_col_names)) param_row[[state_col_names[i]]] <- est_betas[i]
            } else {
               for(n in names(est_betas)) param_row[[paste0("state_", n)]] <- est_betas[[n]]
            }
            # 2. Obs
            if(length(est_alphas) == length(obs_col_names)){
               for(i in seq_along(obs_col_names)) param_row[[obs_col_names[i]]] <- est_alphas[i]
            } else {
               for(n in names(est_alphas)) param_row[[paste0("obs_", n)]] <- est_alphas[[n]]
            }
            
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
final_stats  <- dplyr::bind_rows(all_dataset_stats)

write.csv(final_params, file.path(main_output_dir, "gradient_parameters.csv"), row.names = FALSE)
write.csv(final_preds, file.path(main_output_dir, "gradient_predictions.csv"), row.names = FALSE)
write.csv(final_stats, file.path(main_output_dir, "dataset_descriptive_stats.csv"), row.names = FALSE)

cat("\n--- Gradient Experiment Complete ---\n")