library(dplyr)
library(rje) # For expit()
library(terra) # Added for raster operations
library(sf) # Added for spatial operations
library(Matrix) # Added for sparse 'w' matrix

source(file.path("R","utils.R"))
source(file.path("R","clustering_helpers.R"))
source(file.path("R","model_helpers.R"))


prepare_train_data <- function (
    state_covs, 
    obs_covs,
    cov_tif,  
    placeholder_spec_name = "AMCR"
){

  train_filename <- paste0(placeholder_spec_name, "_zf_filtered_region_2017.csv")
  train_df_og <- read.delim(
    file.path("checklist_data","species", placeholder_spec_name, train_filename), 
    sep = ",", header = T
  )

  train_df_og <- train_df_og[!is.na(train_df_og$duration_minutes),]
  train_df_og <- train_df_og[
    train_df_og$observation_date >= "2017-05-15" & 
    train_df_og$observation_date <= "2017-07-09",
  ]
   
  train_df <- train_df_og

  # Use the standardized function name from utils.R
  train_env_df <- extract_state_covs(train_df, cov_tif) 

  # Now, join the original data with the new, numeric covariate data
  train_df <- inner_join(train_df, train_env_df, by = "checklist_id")
  
  # Use standardized var names
  norm_res <- norm_ds(train_df, obs_covs, state_covs) 

  train_df <- norm_res$df
  norm_list <- norm_res$n_l

  train_df$species_observed <- -1
  train_df$occupied_prob <- -1
  train_df$det_prob <- -1
   
  train_df$formatted_date <- train_df$observation_date

  return (list(train_df = train_df, norm_list = norm_list))
}


# +++ NEW FUNCTION +++
prepare_test_data <- function (
    state_covs, 
    obs_covs,
    cov_tif,
    norm_list, # Pass in the norm_list from the training data
    placeholder_spec_name = "AMCR"
){
  
  test_filename <- paste0(placeholder_spec_name, "_zf_filtered_region_2018.csv")
  test_df_og <- read.delim(
    file.path("checklist_data","species", placeholder_spec_name, test_filename), 
    sep = ",", header = T
  )
  
  test_df_og <- test_df_og[!is.na(test_df_og$duration_minutes),]
  test_df_og <- test_df_og[
    test_df_og$observation_date >= "2018-05-15" & 
    test_df_og$observation_date <= "2018-07-09",
  ]
  test_df <- test_df_og
  
  # Use the standardized function name from utils.R
  test_env_df <- extract_state_covs(test_df, cov_tif) 
  
  test_df <- inner_join(test_df, test_env_df, by = "checklist_id")
  
  # Use standardized var names AND the norm_list from training
  norm_res <- norm_ds(test_df, obs_covs, state_covs, norm_list = norm_list) 
  
  test_df <- norm_res$df
  
  return (test_df) # Just return the processed dataframe
}


# +++ MODIFIED FUNCTION +++
simulate_train_data <-  function (
    reference_clustering_df,
    site_geoms_sf,
    # parameter_set_row,  # Removed
    # state_cov_names,    # Removed
    obs_cov_names,
    obs_par_list,         # Added
    N_j_raster            # Added
) {
  
  # === 1. EXTRACT PARAMETERS ===
  message("  (sim_train) Extracting parameters...")
  # state_par_list no longer needed
  # obs_par_list is now passed directly
  
  # === 2. PROJECT RASTER... (REMOVED) ===
  
  # === 3. CALCULATE CELL-LEVEL INTENSITY... (REMOVED) ===
  
  # === 4. CALCULATE CELL-LEVEL EXPECTED ABUNDANCE... (REMOVED) ===
  # N_j_raster is now passed directly
  
  # === 5. CALCULATE SITE-LEVEL EXPECTED ABUNDANCE (lambda_tilde_i) ===
  message("  (sim_train) Extracting site-level expected abundance (lambda_tilde_i)...")
  
  # Use exact=TRUE (implies weights=TRUE) and fun="sum"
  # This sums the cell values (N_j) weighted by the polygon overlap fraction
  lambda_tilde_i_df <- terra::extract(
    N_j_raster,
    site_geoms_sf,
    fun = "sum",
    exact = TRUE,
    ID = FALSE # We already have site IDs in site_geoms_sf
  )
  
  # Store this in the site_geoms_sf dataframe
  site_geoms_sf$lambda_tilde_i <- lambda_tilde_i_df[,1]
  
  # === 6. SIMULATE SITE-LEVEL ABUNDANCE (N_i) & OCCUPANCY (Z_i) ===
  message("  (sim_train) Simulating site-level abundance (N_i) and occupancy (Z_i)...")
  
  # \psi_i = 1 - e^{-\tilde{\lambda}_i} 
  site_geoms_sf$psi_i <- 1 - exp(-site_geoms_sf$lambda_tilde_i)
  
  # N_i ~ Poisson(\tilde{\lambda}_i) 
  site_geoms_sf$N_i <- rpois(
    n = nrow(site_geoms_sf),
    lambda = site_geoms_sf$lambda_tilde_i
  )
  
  # Z_i = 1 if N_i > 0, 0 otherwise
  site_geoms_sf$Z_i <- ifelse(site_geoms_sf$N_i > 0, 1, 0)
  
  
  # === 7. MERGE SITE-LEVEL STATE TO CHECKLISTS ===
  message("  (sim_train) Merging site state to checklists...")
  
  # Create a lookup table from the sf object
  # *** ADDED N_i TO THE LOOKUP ***
  site_state_lookup <- sf::st_drop_geometry(
    site_geoms_sf[, c("site", "psi_i", "Z_i", "N_i")]
  )
  names(site_state_lookup) <- c("site", "occupied_prob", "occupied", "N")
  
  # Join with the original checklist data
  res_df <- dplyr::left_join(
    reference_clustering_df, 
    site_state_lookup, 
    by = "site"
  )
  
  # === 8. SIMULATE CHECKLIST-LEVEL DETECTION (y_it) ===
  # y_it | Z_i ~ Bernoulli(Z_i * p_it) 
  message("  (sim_train) Simulating checklist-level detection...")
  
  # Calculate detection probability p_it for all checklists
  det_logit <- calculate_weighted_sum(obs_par_list, res_df)
  res_df$det_prob <- rje::expit(det_logit)
  
  # Simulate a detection attempt for every checklist
  res_df$detection <- rbinom(
    n = nrow(res_df), 
    size = 1, 
    prob = res_df$det_prob
  )
  
  # === 9. FINALIZE OBSERVATION ===
  # species_observed = Z_i * detection
  res_df$species_observed <- res_df$occupied * res_df$detection
  
  message("  (sim_train) Simulation of training data complete.")
  
  # Keep columns consistent with old code and test data simulation
  # (e.g., `occupied`, `det_prob`, `species_observed`, `N`)
  return (res_df)
}


# +++ CORRECTED FUNCTION +++
simulate_test_data <- function (
    base_test_df,
    obs_cov_names,
    obs_par_list,
    N_j_raster,
    albers_crs_str,
    area_j_raster  # <--- ADD THIS ARGUMENT
){
  
  # === 1. START WITH PRE-PROCESSED DATA ===
  test_df <- base_test_df
  
  # === 2. EXTRACT PARAMETERS (Observation only) ===
  # state_par_list is no longer needed
  
  # === 3. PROJECT TEST POINTS ===
  test_sf <- sf::st_as_sf(
    test_df, 
    coords = c("longitude", "latitude"), 
    crs = "+proj=longlat +datum=WGS84"
  )
  test_sf_albers <- sf::st_transform(test_sf, crs = albers_crs_str)
  test_df$area_j <- terra::extract(area_j_raster, test_sf_albers, ID = FALSE)[,1]
  
  # === 4. SIMULATE ABUNDANCE (N_j) & OCCUPANCY (Z_j) (CELL-LEVEL) ===
  
  # +++ FIX: Extract expected abundance directly from N_j_raster +++
  # This is the expected abundance for the cell
  lambda_j_cell_df <- terra::extract(N_j_raster, test_sf_albers, ID = FALSE)
  test_df$lambda_j_cell <- lambda_j_cell_df[, 1]

  # Simulate N_j ~ Poisson(lambda_j_cell)
  test_df$N <- rpois(n = nrow(test_df), lambda = test_df$lambda_j_cell)
  
  # Derive occupied state Z_j from N_j
  test_df$Z_i <- ifelse(test_df$N > 0, 1, 0)
  
  # Store the true occupancy probability (P(N_j > 0))
  test_df$occupied_prob <- 1 - exp(-test_df$lambda_j_cell)

  # === 5. SIMULATE DETECTION (y_j) ===
  test_df$det_prob <- calculate_weighted_sum(obs_par_list, test_df)
  test_df$det_prob <- rje::expit(test_df$det_prob)
  test_df$detection <- rbinom(nrow(test_df), 1, test_df$det_prob)
  
  # === 6. FINAL OBSERVATION ===
  test_df$species_observed <- test_df$Z_i * test_df$detection

  message("  (sim_test) Simulation of test data complete.")

  return (test_df)
}