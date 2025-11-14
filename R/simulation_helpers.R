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


simulate_train_data <-  function (
    reference_clustering_df,
    site_geoms_sf,
    parameter_set_row, 
    state_cov_names, 
    obs_cov_names,
    cov_tif,
    norm_list
) {
  
  # === 1. EXTRACT PARAMETERS ===
  message("  (sim_train) Extracting parameters...")
  state_par_list <- as.list(parameter_set_row[, c("state_intercept", state_cov_names)])
  obs_par_list <- as.list(parameter_set_row[, c("obs_intercept", obs_cov_names)])
  
  names(state_par_list)[1] <- "intercept"
  names(obs_par_list)[1] <- "intercept"
  
  # === 2. PROJECT RASTER TO MATCH SITE GEOMETRIES ===
  message("  (sim_train) Projecting covariate raster...")
  # Get the Albers equal-area CRS from the site geometries
  albers_crs_str <- sf::st_crs(site_geoms_sf)$wkt
  
  # Project the covariate raster to this CRS
  # This ensures area calculations are in consistent units (meters)
  cov_tif_albers <- terra::project(cov_tif, albers_crs_str, method="bilinear")
  
  # === 3. CALCULATE CELL-LEVEL INTENSITY (lambda_j) ===
  # This implements \lambda(s) = e^{f(x_s)} 
  message("  (sim_train) Calculating cell-level intensity raster (lambda_j)...")
  
  # Start with the intercept
  log_lambda_j_raster <- cov_tif_albers[[1]] * 0 + state_par_list$intercept
  
  # Add weighted covariates
  for (cov_name in state_cov_names) {
    if (cov_name %in% names(cov_tif_albers)) {
      log_lambda_j_raster <- log_lambda_j_raster + (cov_tif_albers[[cov_name]] * state_par_list[[cov_name]])
    }
  }
  
  # lambda_j = exp(f(x_j))
  lambda_j_raster <- exp(log_lambda_j_raster)
  
  # === 4. CALCULATE CELL-LEVEL EXPECTED ABUNDANCE (N_j) ===
  # N_j = \lambda_j * Area_j
  message("  (sim_train) Calculating cell-level expected abundance (N_j)...")
  
  # Get cell area in square meters
  area_j_raster <- terra::cellSize(cov_tif_albers, unit="m")
  
  # N_j_raster now holds the expected number of individuals per cell
  N_j_raster <- lambda_j_raster * area_j_raster
  
  # === 5. CALCULATE SITE-LEVEL EXPECTED ABUNDANCE (lambda_tilde_i) ===
  # This implements \tilde{\lambda}_i \approx \sum_j \lambda_j \cdot Area(B_i \cap C_j) 
  # by extracting from the N_j raster: \sum_j (N_j * fraction_ij)
  # which is \sum_j ( (\lambda_j * Area_j) * (Area_ij / Area_j) ) = \sum_j (\lambda_j * Area_ij)
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
  
  # === 6. SIMULATE SITE-LEVEL OCCUPANCY (Z_i) ===
  # \psi_i = 1 - e^{-\tilde{\lambda}_i} 
  message("  (sim_train) Simulating site-level occupancy (Z_i)...")
  
  site_geoms_sf$psi_i <- 1 - exp(-site_geoms_sf$lambda_tilde_i)
  
  # Z_i ~ Bernoulli(\psi_i) 
  site_geoms_sf$Z_i <- rbinom(
    n = nrow(site_geoms_sf),
    size = 1,
    prob = site_geoms_sf$psi_i
  )
  
  # === 7. MERGE SITE-LEVEL STATE TO CHECKLISTS ===
  message("  (sim_train) Merging site state to checklists...")
  
  # Create a lookup table from the sf object
  site_state_lookup <- sf::st_drop_geometry(
    site_geoms_sf[, c("site", "psi_i", "Z_i")]
  )
  names(site_state_lookup) <- c("site", "occupied_prob", "occupied")
  
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
  # (e.g., `occupied`, `det_prob`, `species_observed`)
  return (res_df)
}



simulate_test_data <- function (
    norm_list, 
    parameter_set_row, 
    state_cov_names, 
    obs_cov_names, 
    cov_tif, # <-- Pass the raster in
    placeholder_spec_name = "AMCR"
){
  
  # === 1. LOAD & FILTER TEST DATA ===

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
  
  # === 2. EXTRACT & NORMALIZE COVARIATES ===
  test_env_df <- extract_state_covs(test_df, cov_tif) # Use new function name
  test_df <- inner_join(test_df, test_env_df, by = "checklist_id")
  
  # Use standardized var names
  norm_res <- norm_ds(test_df, obs_cov_names, state_cov_names, norm_list = norm_list)
  
  test_df <- norm_res$df
  
  # === 3. EXTRACT PARAMETERS ===
  state_par_list <- as.list(parameter_set_row[, c("state_intercept", state_cov_names)])
  obs_par_list <- as.list(parameter_set_row[, c("obs_intercept", obs_cov_names)])
  
  names(state_par_list)[1] <- "intercept"
  names(obs_par_list)[1] <- "intercept"
  

  # === 4. SIMULATE OCCUPANCY & DETECTION (CHECKLIST-LEVEL) ===
  # For test data, we simulate everything at the checklist-level
  
  test_df$occupied_prob <- calculate_weighted_sum(state_par_list, test_df)
  test_df$occupied_prob <- rje::expit(test_df$occupied_prob)
  test_df$occupied <- rbinom(nrow(test_df), 1, test_df$occupied_prob)

  test_df$det_prob <- calculate_weighted_sum(obs_par_list, test_df)
  test_df$det_prob <- rje::expit(test_df$det_prob)
  test_df$detection <- rbinom(nrow(test_df), 1, test_df$det_prob)

  
  # === 5. FINAL OBSERVATION ===
  test_df$species_observed <- test_df$occupied * test_df$detection

  message("  (sim_test) Simulation of test data complete.")

  return (test_df)
}