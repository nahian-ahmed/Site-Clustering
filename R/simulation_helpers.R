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
    site_geoms_sf,
    parameter_set_row, 
    state_cov_names, 
    obs_cov_names,
    cov_tif,
    norm_list
) {
  
  # === 1. EXTRACT PARAMETERS ===
  cat("    - (Simulating occuN) Extracting parameters...\n")
  state_par_list <- as.list(parameter_set_row[, c("state_intercept", state_cov_names)])
  obs_par_list <- as.list(parameter_set_row[, c("obs_intercept", obs_cov_names)])
  
  names(state_par_list)[1] <- "intercept"
  names(obs_par_list)[1] <- "intercept"
  

  # === 2. GET CELL-LEVEL COVARIATES (cellCovs) ===
  cat("    - (Simulating occuN) Extracting cellCovs from raster...\n")
  cell_data <- terra::as.data.frame(cov_tif, xy = TRUE, cells = TRUE)
  colnames(cell_data)[colnames(cell_data) == 'cell'] <- 'cell_id'
  
  norm_res <- norm_ds(cell_data, obs_covs = NULL, state_covs = state_cov_names, norm_list = norm_list)
  cellCovs <- norm_res$df
  n_cells <- nrow(cellCovs)


  # === 3. CALCULATE CELL-LEVEL LAMBDA (lambda_j) ===
  cat("    - (Simulating occuN) Calculating cell-level lambda...\n")
  cellCovs$lambda_j <- exp(calculate_weighted_sum(state_par_list, cellCovs))


  # === 4. DEFINE SITES & BUILD 'w' MATRIX (NEW AREAL OVERLAP LOGIC) ===
  
  albers_crs <- sf::st_crs(site_geoms_sf)
  
  # 2. Get raster CRS and calculate cell areas in m^2
  cat("    - (Simulating occuN) Calculating cell areas...\n")
  rast_crs <- sf::st_crs(cov_tif)
  # Calculate the area of each cell in m^2
  cell_areas_m2_rast <- terra::expanse(cov_tif, unit="m", transform=TRUE)
  cell_areas_df <- terra::as.data.frame(cell_areas_m2_rast, cells=TRUE)
  # Rename 'area' column to avoid conflicts
  colnames(cell_areas_df)[colnames(cell_areas_df) == 'area'] <- 'cell_area_m2'

  # 3. Transform site polygons to match the raster's CRS for extraction
  cat("    - (Simulating occuN) Transforming site geometries to raster CRS...\n")
  site_geoms_wgs84 <- sf::st_transform(site_geoms_sf, crs = rast_crs)
  
  # 4. Extract areal overlap weights using terra::extract
  # This is the memory-safe replacement for st_intersection
  cat("    - (Simulating occuN) Extracting areal overlap weights (using terra::extract)...\n")
  site_cell_map_raw <- terra::extract(
    cov_tif,
    site_geoms_wgs84,
    exact = TRUE,      # <-- This calculates fractional overlap
    weights = TRUE,    # <-- This returns the weights
    ID = TRUE         # We will use our own site IDs
  )
  
  # 5. Map results back to site and cell IDs
  # 'ID' from extract() is the index of the polygon (1 to 625)
  # 'cell' is the cell ID
  # 'weight' is the fractional overlap (0 to 1)
  site_cell_map_raw$site <- site_geoms_wgs84$site[site_cell_map_raw$ID]
  
  # 6. Calculate the final overlap area in m^2 (fraction * total_cell_area)
  cat("    - (Simulating occuN) Calculating overlap areas in m^2...\n")
  site_cell_map_with_areas <- merge(site_cell_map_raw, cell_areas_df, by.x="cell", by.y="cell")
  site_cell_map_with_areas$overlap_area_m2 <- site_cell_map_with_areas$weight * site_cell_map_with_areas$cell_area_m2

  # 7. Create the final map for the sparse matrix
  site_cell_map <- site_cell_map_with_areas[, c("site", "cell", "overlap_area_m2")]
  colnames(site_cell_map)[colnames(site_cell_map) == 'cell'] <- 'cell_id'
  
  # Make sure site IDs are consistent
  unique_site_ids <- unique(site_geoms_sf$site)
  M <- length(unique_site_ids)
  
  cat(sprintf("    - (Simulating occuN) Building %d x %d 'w' matrix with area weights...\n", M, n_cells))
  
  w <- Matrix::sparseMatrix(
      i = match(site_cell_map$site, unique_site_ids),
      j = match(site_cell_map$cell_id, cellCovs$cell_id),
      x = site_cell_map$overlap_area_m2, # <-- THIS IS THE KEY CHANGE
      dims = c(M, n_cells),
      dimnames = list(unique_site_ids, cellCovs$cell_id)
  )
  
  # --- IMPORTANT NOTE ON UNITS ---
  # Your w_ij values are now in m^2[cite: 121].
  # This means the latent abundance (lambda_j) is interpreted as
  # "expected individuals per m^2" in that cell.
  # And lambda_tilde_i is the total expected individuals for the site.
  # This is all 100% correct according to the theory.

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
  cat("    - (Simulating occuN) Enforcing closure and simulating detection...\n")
  # This should be reference_clustering_df, not sites_df
  sites_list <- unique(reference_clustering_df$site)
  closed_df <- enforceClosure(reference_clustering_df, state_cov_names, sites_list)
  
  
  # === 7. MAP SITE-LEVEL OCCUPANCY TO CHECKLISTS ===
  res_df <- dplyr::left_join(closed_df, site_lookup, by = "site")
  
  
  # === 8. SIMULATE DETECTION (CHECKLIST-LEVEL) ===
  res_df$det_prob <- calculate_weighted_sum(obs_par_list, res_df)
  res_df$det_prob <- rje::expit(res_df$det_prob)
  res_df$detection <- rbinom(nrow(res_df), 1, res_df$det_prob)
  
  
  # === 9. FINAL OBSERVATION ===
  res_df$species_observed <- res_df$occupied * res_df$detection
  
  # Add a check for NAs in occupancy, which can happen if a site
  # somehow didn't overlap any cells (e.g., in a lake/NA area)
  n_na_occ <- sum(is.na(res_df$occupied))
  if (n_na_occ > 0) {
      cat(sprintf("    - (Simulating occuN) WARNING: %d checklists had NA occupancy (no cell overlap). Setting to 0.\n", n_na_occ))
      res_df$occupied[is.na(res_df$occupied)] <- 0
      res_df$species_observed[is.na(res_df$species_observed)] <- 0
  }
  
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