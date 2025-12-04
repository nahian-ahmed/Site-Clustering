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
  scale_res <- standardize_ds(train_df, obs_covs, state_covs = NULL) 

  train_df <- scale_res$df
  standardization_params <- scale_res$standardization_params

  train_df$species_observed <- -1
  train_df$occupied_prob <- -1
  train_df$det_prob <- -1
   
  train_df$formatted_date <- train_df$observation_date

  return (list(train_df = train_df, standardization_params = standardization_params))
}


# +++ NEW FUNCTION +++
prepare_test_data <- function (
    state_covs, 
    obs_covs,
    cov_tif,
    standardization_params, # Pass in the norm_list from the training data
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
  scale_res <- standardize_ds(test_df, obs_covs, state_covs = NULL, standardization_params =  standardization_params) 
  
  test_df <- scale_res$df
  
  return (test_df) # Just return the processed dataframe
}





# R/simulation_helpers.R

simulate_train_data <- function (
    reference_clustering_df,
    obs_cov_names,
    obs_par_list,
    w_matrix,             # <-- NEW: Pass the sparse matrix (Sites x Cells)
    cell_density_vector  # <-- NEW: Pass Density values (exp(log_lambda)), NOT Abundance (N_j)
) {
  
  # === 1. CALCULATE SITE-LEVEL EXPECTED ABUNDANCE (lambda_tilde_i) ===
  # Matrix Multiplication: [Sites x Cells] * [Cells x 1] = [Sites x 1]
  # w_matrix contains 'Area overlap', so: Area * Density = Abundance
  
  # Ensure NAs in density are 0 to prevent propagation errors
  cell_density_vector[is.na(cell_density_vector)] <- 0
  
  lambda_tilde_i <- as.numeric(w_matrix %*% cell_density_vector)
  
  # Create a lookup for site state
  # We assume rownames(w_matrix) are the site IDs
  site_state_df <- data.frame(
    site = rownames(w_matrix),
    lambda_tilde_i = lambda_tilde_i
  )
  
  # === 2. SIMULATE SITE-LEVEL ABUNDANCE (N_i) & OCCUPANCY (Z_i) ===
  # psi_i = 1 - exp(-lambda)
  site_state_df$occupied_prob <- 1 - exp(-site_state_df$lambda_tilde_i)
  
  # N_i ~ Poisson(lambda)
  site_state_df$N <- rpois(n = nrow(site_state_df), lambda = site_state_df$lambda_tilde_i)
  
  # Z_i = 1 if N > 0
  site_state_df$occupied <- ifelse(site_state_df$N > 0, 1, 0)
  
  
  # === 3. MERGE TO CHECKLISTS ===
  reference_clustering_df$site <- as.character(reference_clustering_df$site)
  site_state_df$site <- as.character(site_state_df$site)
  
  res_df <- dplyr::left_join(
    reference_clustering_df, 
    site_state_df, 
    by = "site"
  )
  
  # === 4. SIMULATE DETECTION (y_it) ===
  # Calculate detection probability p_it
  det_logit <- calculate_weighted_sum(obs_par_list, res_df)
  res_df$det_prob <- rje::expit(det_logit)
  
  # Simulate detection
  res_df$detection <- rbinom(
    n = nrow(res_df), 
    size = 1, 
    prob = res_df$det_prob
  )
  
  # === 5. FINALIZE ===
  res_df$species_observed <- res_df$occupied * res_df$detection
  
  return (res_df)
}


simulate_test_data <- function (
    base_test_df,
    obs_cov_names,
    obs_par_list,
    w_matrix,             # <-- NEW: Sparse Matrix
    cell_density_vector   # <-- NEW: Density (exp(log_lambda)), NOT N_j_raster
){
  
  test_df <- base_test_df
  
  # === 1. Calculate Site-Level Expected Abundance (Lambda) ===
  # Matrix Mult: [TestSites x Cells] * [Cells x 1] = [TestSites x 1]
  # W matrix contains Area (km^2). Density is per km^2. Result is Abundance.
  
  cell_density_vector[is.na(cell_density_vector)] <- 0
  lambda_vals <- as.numeric(w_matrix %*% cell_density_vector)
  
  test_df$lambda_j_site <- lambda_vals
  
  # === 2. Simulate Biological State ===
  test_df$N <- rpois(n = nrow(test_df), lambda = test_df$lambda_j_site)
  test_df$Z_i <- ifelse(test_df$N > 0, 1, 0)
  test_df$occupied_prob <- 1 - exp(-test_df$lambda_j_site)

  # === 3. Simulate Detection ===
  test_df$det_prob <- calculate_weighted_sum(obs_par_list, test_df)
  test_df$det_prob <- rje::expit(test_df$det_prob)
  test_df$detection <- rbinom(nrow(test_df), 1, test_df$det_prob)
  
  # === 4. Finalize ===
  test_df$species_observed <- test_df$Z_i * test_df$detection

  # Ensure standard dataframe (should be already, but just in case)
  if (inherits(test_df, "sf")) {
    test_df <- sf::st_drop_geometry(test_df)
  }

  return (test_df)
}


# R/simulation_helpers.R

#' Generate Spatial Structures for Test Data
#' 
#' Creates Voronoi geometries, calculates site areas, and builds the sparse W matrix
#' for the test dataset.
#'
#' @param test_df Dataframe containing test coordinates (must have 'longitude', 'latitude').
#' @param albers_crs String. The projection CRS (e.g., Albers) to use for geometry generation.
#' @param buffer_m Numeric. Buffer distance in meters.
#' @param cov_raster_albers SpatRaster. The reference raster projected to albers_crs.
#' @param area_raster SpatRaster. The cell size raster (usually in km^2 or m^2).
#'
#' @return A list containing:
#'   \item{test_df}{The updated dataframe with 'site' IDs and 'area_j' columns.}
#'   \item{w_matrix}{The sparse matrix (Sites x Cells) representing spatial overlap.}
prepare_test_spatial_structures <- function(test_df, 
                                            albers_crs, 
                                            buffer_m, 
                                            cov_raster_albers, 
                                            area_raster) {
  
  cat("--- Generating test geometries and W matrix... ---\n")
  
  # 1. Assign Site IDs
  test_df$site <- seq_len(nrow(test_df))
  
  # 2. Create SF object
  test_sf <- sf::st_as_sf(
    test_df, 
    coords = c("longitude", "latitude"), 
    crs = "+proj=longlat +datum=WGS84"
  )
  
  # 3. Transform and Buffer
  test_sf_albers <- sf::st_transform(test_sf, crs = albers_crs)
  
  # Note: Ensure voronoi_clipped_buffers is available (sourced from model_helpers.R)
  test_geoms <- voronoi_clipped_buffers(test_sf_albers, buffer_dist = buffer_m)
  
  # 4. Calculate Site Area (Static)
  # Sum of raster cell areas within the buffer
  # exact=TRUE calculates fraction of cell covered
  test_area_vals <- terra::extract(area_raster, test_geoms, fun = "sum", exact = TRUE, ID = FALSE)
  test_df$area_j <- test_area_vals[, 1]
  
  # 5. Create W Matrix
  
  # Convert to SpatVector for extraction
  test_vect <- terra::vect(test_geoms)
  test_vect_proj <- terra::project(test_vect, terra::crs(cov_raster_albers))
  
  # Extract overlap (cells, exact fraction, and ID)
  test_overlap_df <- terra::extract(
    cov_raster_albers[[1]], 
    test_vect_proj, 
    cells = TRUE, 
    exact = TRUE, 
    ID = TRUE
  )
  
  # Calculate Area Weights
  test_overlap_cell_areas <- terra::extract(area_raster, test_overlap_df$cell)
  # test_overlap_df$w_area <- test_overlap_df$fraction * test_overlap_cell_areas[,1]
  test_overlap_df$w_area <- test_overlap_df$fraction
  
  # Construct Sparse Matrix (Rows = Test Sites, Cols = Raster Cells)
  n_test_sites <- nrow(test_geoms)
  n_cells <- terra::ncell(cov_raster_albers)
  
  w_matrix <- Matrix::sparseMatrix(
    i = test_overlap_df$ID,
    j = test_overlap_df$cell,
    x = test_overlap_df$w_area,
    dims = c(n_test_sites, n_cells)
  )
  
  cat("--- Test spatial structures created. ---\n")
  
  # Return updated DF and Matrix
  return(list(
    test_df = test_df,
    w_matrix = w_matrix
  ))
}