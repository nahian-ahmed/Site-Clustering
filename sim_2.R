# -----------------------------------------------------------------
# sim_2.R
# Complexity Gradient: Bridging Simulation and Reality
# 4 Variants: Uniform/Real Locs x Simple/Complex Covariates
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

# Source existing helpers
source(file.path("R", "utils.R"))
source(file.path("R", "simulation_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))
# Note: Ensure plotting_helpers.R is sourced if you use plotting functions
# source(file.path("R", "plotting_helpers.R")) 

set.seed(123)

###
# 2. CONFIGS & VARIANTS
###

# Load Parameters and Clusterings from CSV (Same as sim_3.R)
sim_params <- read.csv(file.path("config", "simulation_parameters.csv"))
sim_clusterings <- read.csv(file.path("config", "simulation_clusterings.csv"))

# Define the 4 Variants
# We Map the CSV column names to our Synthetic Columns here
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
n_simulations <- 10  # Adjust as needed
n_fit_repeats <- 10
selected_optimizer <- "nlminb"
buffer_m <- 200
res_m <- 100

# Output Directory
output_dir <- file.path("simulation_experiments", "output", "sim_2")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

###
# 3. PREPROCESS RASTER (Global)
###

# Standardize Raster Once (Same as sim_3.R)
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = res_m)

# Standardize
standardization_results <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- standardization_results$raster
state_cov_params <- standardization_results$params

# Area Raster
cell_area_km2 <- (res_m / 1000) * (res_m / 1000)
area_j_raster <- cov_tif_albers[[1]] * 0 + cell_area_km2

# Full Covariate Dataframe (for fitting)
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))
full_raster_covs[is.na(full_raster_covs)] <- 0


###
# 4. HELPER: GENERATE DATASET FOR VARIANT
###
generate_variant_dataset <- function(variant_cfg, cov_raster, n_target_obs, real_file_template) {
  
  # A. Locations
  if (variant_cfg$loc_type == "uniform") {
    cat("    [DataGen] Sampling Uniform Locations...\n")
    # Sample from raster cells (ensuring valid data)
    valid_cells <- terra::cells(cov_raster[[1]])
    # Sample indices
    sampled_indices <- sample(valid_cells, n_target_obs, replace = TRUE) # Replace=TRUE allows coincident points like eBird
    
    coords_proj <- terra::xyFromCell(cov_raster[[1]], sampled_indices)
    
    # Create DF with projected coords, then convert to Lat/Lon for compatibility with helpers
    df_proj <- data.frame(x = coords_proj[,1], y = coords_proj[,2])
    v_proj <- terra::vect(df_proj, geom=c("x", "y"), crs = terra::crs(cov_raster))
    v_geo <- terra::project(v_proj, "+proj=longlat +datum=WGS84")
    coords_geo <- terra::crds(v_geo)
    
    base_df <- data.frame(
      checklist_id = paste0("unif_", 1:n_target_obs),
      locality_id = paste0("loc_unif_", 1:n_target_obs), # FIX: Added locality_id
      latitude = coords_geo[,2],
      longitude = coords_geo[,1],
      observation_date = "2017-06-01",
      formatted_date = "2017-06-01"
    )
  } else {
    cat("    [DataGen] Using Real eBird Locations...\n")
    # Load Real Data (using sim_3.R logic)
    real_df <- read.csv(real_file_template)
    # Filter basic validity
    real_df <- real_df[!is.na(real_df$latitude) & !is.na(real_df$longitude),]
    
    # FIX: Included locality_id in the subset
    base_df <- real_df[, c("checklist_id", "locality_id", "latitude", "longitude", "observation_date")]
    base_df$formatted_date <- base_df$observation_date
    
    # Ensure N matches if needed, but usually we just take the full filtered set
    if(nrow(base_df) > n_target_obs) base_df <- base_df[1:n_target_obs, ]
  }
  
  # B. State Covariates (Extract from Raster)
  # We extract ALL state covs for the dataframe, but will only USE specific ones in Sim
  env_df <- extract_state_covs(base_df, cov_raster)
  base_df <- dplyr::inner_join(base_df, env_df, by = "checklist_id")
  
  # C. Observation Covariates (Simulated)
  # We create synthetic columns: obs_cov_1, obs_cov_2...
  cat(sprintf("    [DataGen] Simulating %d Obs Covariates...\n", variant_cfg$n_obs_covs))
  
  synthetic_obs_names <- c()
  for(i in 1:variant_cfg$n_obs_covs){
    col_name <- paste0("obs_sim_", i)
    base_df[[col_name]] <- rnorm(nrow(base_df), mean = 0, sd = 1)
    synthetic_obs_names <- c(synthetic_obs_names, col_name)
  }
  
  return(list(df = base_df, obs_cols = synthetic_obs_names))
}


###
# 5. MAIN LOOP: VARIANTS
###

all_param_results <- list()

# Get target N from a real file (to match eBird N)
template_file <- file.path("checklist_data", "species", "AMCR", "AMCR_zf_filtered_region_2017.csv")
temp_df <- read.csv(template_file)
TARGET_N_OBS <- nrow(temp_df)


for (v_name in names(variants)) {
  
  variant <- variants[[v_name]]
  cat(paste("\n################################################\n"))
  cat(paste("### STARTING VARIANT:", v_name, "###\n"))
  cat(paste("################################################\n"))
  
  # 1. Generate Dataset for this Variant (Fixed Geometry for the Variant)
  #    We generate ONE dataset and hold it constant across simulations (like a real dataset)
  ds_bundle <- generate_variant_dataset(variant, cov_tif_albers, TARGET_N_OBS, template_file)
  current_train_df <- ds_bundle$df
  current_obs_cols <- ds_bundle$obs_cols
  current_state_cols <- variant$state_covs_used
  
  # 2. Pre-compute Clusterings (Using sim_3.R logic)
  #    We only need the Reference Clustering defined in sim_clusterings.csv
  #    But we calculate them all for consistency if comparison is needed
  cat("  --- Pre-computing Clusterings for Variant ---\n")
  # Note: get_clusterings expects 'og_data' and 'state_covs'
  # We pass ALL state covs for clustering to work, even if Sim uses fewer
  all_clusterings <- get_clusterings(unique(sim_clusterings$method), current_train_df, names(full_raster_covs))
  
  # 3. Pre-compute Site Geometries & W Matrices
  cat("  --- Generating Geometries & W Matrices ---\n")
  all_w_matrices <- list()
  
  for (m_name in names(all_clusterings)) {
    cluster_data <- all_clusterings[[m_name]]
    if (is.list(cluster_data)) cluster_data <- cluster_data$result_df
    
    # Geometries
    geoms <- create_site_geometries(cluster_data, cov_tif_albers, buffer_m, m_name)
    
    # Disjoint Split (Important!)
    split_res <- disjoint_site_geometries(geoms, cluster_data)
    geoms <- split_res$geoms
    all_clusterings[[m_name]] <- split_res$data # Update DF
    
    # W Matrix
    all_w_matrices[[m_name]] <- generate_overlap_matrix(geoms, cov_tif_albers)
  }
  
  
  # 4. Simulation Loops (Params -> Sims)
  for (param_idx in seq_len(nrow(sim_params))) {
    
    # --- A. Setup Parameters ---
    curr_params <- sim_params[param_idx, ]
    
    # Filter State Params: Force unused state covs to 0 for the "Truth"
    state_par_list <- list(intercept = curr_params$state_intercept)
    for(sc in names(full_raster_covs)){
      if(sc %in% current_state_cols){
        state_par_list[[sc]] <- curr_params[[sc]]
      } else {
        state_par_list[[sc]] <- 0 # Zero out unused covs
      }
    }
    
    # Map Obs Params: Map the first N params from CSV to our N synthetic columns
    # CSV headers: duration_minutes, effort_distance_km, etc.
    obs_csv_headers <- c("duration_minutes", "effort_distance_km", "number_observers", "time_observations_started", "day_of_year")
    
    obs_par_list <- list(intercept = curr_params$obs_intercept)
    for(i in 1:variant$n_obs_covs){
      # Take coefficient from the i-th CSV column
      csv_col <- obs_csv_headers[i]
      # Assign to the i-th Synthetic column
      obs_par_list[[ current_obs_cols[i] ]] <- curr_params[[csv_col]]
    }
    
    # --- B. Calculate Density Surface (Truth) ---
    log_lambda_j <- cov_tif_albers[[1]] * 0 + state_par_list$intercept
    for (nm in names(full_raster_covs)) {
        if(state_par_list[[nm]] != 0){
             log_lambda_j <- log_lambda_j + (cov_tif_albers[[nm]] * state_par_list[[nm]])
        }
    }
    cell_density_val <- terra::values(exp(log_lambda_j), mat = FALSE)
    cell_density_val[is.na(cell_density_val)] <- 0
    
    
    # --- C. Loop Simulations ---
    ref_method <- sim_clusterings$method[1] # Assuming 1 ref method for simplicity, or loop rows
    
    for(sim_num in 1:n_simulations){
      set.seed(sim_num * param_idx) # Unique seed
      
      cat(sprintf("  [V: %s | P: %d | Sim: %d] Ref: %s\n", v_name, param_idx, sim_num, ref_method))
      
      # Simulate Data
      # We use the Reference Clustering DF and W Matrix
      ref_df <- all_clusterings[[ref_method]]
      if (is.list(ref_df)) ref_df <- ref_df$result_df
      ref_w <- all_w_matrices[[ref_method]]
      
      train_data <- simulate_train_data(
        reference_clustering_df = ref_df,
        obs_cov_names = current_obs_cols, # Pass synthetic names
        obs_par_list = obs_par_list,
        w_matrix = ref_w,
        cell_density_vector = cell_density_val
      )
      
      # Fit Model (Using Reference Method for now, or loop 'comparison_method_list')
      # We just fit the Reference Method to check recovery
      
      umf <- prepare_occuN_data(train_data, ref_df, ref_w, current_obs_cols, full_raster_covs)
      
      # Formulas based on USED covariates
      state_form <- as.formula(paste("~", paste(current_state_cols, collapse = " + ")))
      obs_form <- as.formula(paste("~", paste(current_obs_cols, collapse = " + ")))
      
      fm <- fit_occuN_model(
        umf, 
        state_form, 
        obs_form, 
        n_reps = n_fit_repeats, 
        stable_reps = 3, # Lower for speed in gradient
        optimizer = selected_optimizer,
        lower = -10, upper = 10
      )
      
      # Save Results
      if(!is.null(fm)){
        est_betas <- coef(fm, 'state')
        est_alphas <- coef(fm, 'det')
        
        # We save raw estimates with names to map back later
        res_row <- data.frame(
          Variant = v_name,
          ParamSet = param_idx,
          Sim = sim_num,
          Type = c(rep("State", length(est_betas)), rep("Obs", length(est_alphas))),
          ParamName = c(names(est_betas), names(est_alphas)),
          Estimate = c(as.numeric(est_betas), as.numeric(est_alphas)),
          NLL = fm@negLogLike
        )
        all_param_results[[length(all_param_results)+1]] <- res_row
      }
      
    } # Sim Loop
  } # Param Loop
} # Variant Loop

# Save Final
final_res <- dplyr::bind_rows(all_param_results)
write.csv(final_res, file.path(output_dir, "gradient_results.csv"), row.names = FALSE)
cat("\n--- Gradient Experiment Complete ---\n")