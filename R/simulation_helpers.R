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


# -----------------------------------------------------------------
# --- THIS IS THE UPDATED FUNCTION ---
# -----------------------------------------------------------------
simulate_train_data <-  function (
    reference_clustering_df, 
    parameter_set_row, 
    state_cov_names, 
    obs_cov_names,
    cov_tif,
    norm_list # <-- NEW ARGUMENT
) {
  
  # === 1. EXTRACT PARAMETERS ===
  # (Same as before)
  state_par_list <- as.list(parameter_set_row[, c("state_intercept", state_cov_names)])
  obs_par_list <- as.list(parameter_set_row[, c("obs_intercept", obs_cov_names)])
  
  names(state_par_list)[1] <- "intercept"
  names(obs_par_list)[1] <- "intercept"
  

  # === 2. GET CELL-LEVEL COVARIATES (cellCovs) ===
  # Extract all cell values from the raster
  cat("    - (Simulating occuN) Extracting cellCovs from raster...\n")
  cell_data <- terra::as.data.frame(cov_tif, xy = TRUE, cells = TRUE)
  colnames(cell_data)[colnames(cell_data) == 'cell'] <- 'cell_id'
  
  # Normalize cell-level covariates using the norm_list from training data
  # We assume norm_ds can handle NULL obs_covs
  norm_res <- norm_ds(cell_data, obs_covs = NULL, state_covs = state_cov_names, norm_list = norm_list)
  cellCovs <- norm_res$df
  n_cells <- nrow(cellCovs)


  # === 3. CALCULATE CELL-LEVEL LAMBDA (lambda_j) ===
  # Calculate latent abundance for *every cell* in the landscape
  cat("    - (Simulating occuN) Calculating cell-level lambda...\n")
  cellCovs$lambda_j <- exp(calculate_weighted_sum(state_par_list, cellCovs))


  # === 4. DEFINE SITES & BUILD 'w' MATRIX ===
  # 'w' maps sites (rows) to cells (columns)
  # We define a site's area as the *set of cells* containing its checklists
  
  sites_df <- reference_clustering_df
  unique_site_ids <- unique(sites_df$site)
  M <- length(unique_site_ids)
  
  cat(sprintf("    - (Simulating occuN) Building %d x %d 'w' matrix...\n", M, n_cells))
  
  # Find which cell each checklist belongs to
  sites_sf <- sf::st_as_sf(sites_df, coords = c("longitude", "latitude"), crs = sf::st_crs(cov_tif))
  sites_df$cell_id <- terra::extract(cov_tif, sites_sf, cell = TRUE)$cell
  
  # Get unique (site, cell) pairs
  site_cell_map <- unique(sites_df[, c("site", "cell_id")])
  site_cell_map <- site_cell_map[!is.na(site_cell_map$cell_id), ]
  
  # Build a sparse w matrix
  w <- Matrix::sparseMatrix(
      i = match(site_cell_map$site, unique_site_ids),
      j = match(site_cell_map$cell_id, cellCovs$cell_id),
      x = 1, # Binary weight: 1 if cell is in site, 0 otherwise
      dims = c(M, n_cells),
      dimnames = list(unique_site_ids, cellCovs$cell_id)
  )

  # === 5. SIMULATE SITE-LEVEL OCCUPANCY (Z_i) ===
  # \tilde{\lambda}_i = \sum_j w_{ij} \lambda_j
  cat("    - (Simulating occuN) Aggregating site-level lambda_tilde...\n")
  lambda_tilde_i <- w %*% cellCovs$lambda_j
  
  # \psi_i = 1 - e^{-\tilde{\lambda}_i}
  psi_i <- 1 - exp(-lambda_tilde_i)
  
  # Simulate true occupancy state (Z_i)
  Z_i <- rbinom(M, 1, psi_i)
  
  # Create lookup for site-level values
  site_lookup <- data.frame(
      site = unique_site_ids, 
      occupied_prob = as.numeric(psi_i), 
      occupied = Z_i
  )

  # === 6. ENFORCE CLOSURE ===
  # (Same as before)
  cat("    - (Simulating occuN) Enforcing closure and simulating detection...\n")
  sites_list <- unique(sites_df$site)
  closed_df <- enforceClosure(sites_df, state_cov_names, sites_list)
  
  
  # === 7. MAP SITE-LEVEL OCCUPANCY TO CHECKLISTS ===
  # (Same as before)
  res_df <- dplyr::left_join(closed_df, site_lookup, by = "site")
  
  
  # === 8. SIMULATE DETECTION (CHECKLIST-LEVEL) ===
  # (Same as before)
  res_df$det_prob <- calculate_weighted_sum(obs_par_list, res_df)
  res_df$det_prob <- rje::expit(res_df$det_prob)
  res_df$detection <- rbinom(nrow(res_df), 1, res_df$det_prob)
  
  
  # === 9. FINAL OBSERVATION ===
  # (Same as before)
  res_df$species_observed <- res_df$occupied * res_df$detection
  
  return (res_df)
}
# -----------------------------------------------------------------
# --- END OF UPDATED FUNCTION ---
# -----------------------------------------------------------------



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

  return (test_df)
}