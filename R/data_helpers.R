library(dplyr)
library(terra)
library(sf)
library(Matrix)

# Source necessary dependencies
# extract_state_covs and standardize_ds are in utils.R
# voronoi_clipped_buffers is in model_helpers.R
source(file.path("R", "utils.R"))
source(file.path("R", "model_helpers.R"))


#' Prepare Training Data (Real Data)
#' 
#' Loads the 2017 dataset, extracts raster covariates, and standardizes them.
#' Unlike the simulation version, this PRESERVES the real species_observed values.
#'
prepare_train_data <- function (
    state_covs, 
    obs_covs,
    cov_tif, 
    state_standardization_params = list(), 
    placeholder_spec_name = "AMCR"
){

  # 1. Load Data
  train_filename <- paste0(placeholder_spec_name, "_zf_filtered_region_2017.csv")
  file_path <- file.path("checklist_data", "species", placeholder_spec_name, train_filename)
  
  if (!file.exists(file_path)) {
    stop(paste("Training file not found:", file_path))
  }

  train_df_og <- read.delim(file_path, sep = ",", header = TRUE)

  # 2. Filter
  train_df_og <- train_df_og[!is.na(train_df_og$duration_minutes),]
  train_df_og <- train_df_og[
    train_df_og$observation_date >= "2017-05-15" & 
    train_df_og$observation_date <= "2017-07-09",
  ]
   
  train_df <- train_df_og

  # 3. Extract Raw Raster Values
  # (Assumes cov_tif is the raw raster if standardizing later, or pre-scaled)
  train_env_df <- extract_state_covs(train_df, cov_tif) 

  # 4. Join Covariates to Checklists
  train_df <- inner_join(train_df, train_env_df, by = "checklist_id")
  
  # 5. Standardize Data
  scale_res <- standardize_ds(
      train_df, 
      obs_covs, 
      state_covs = state_covs, 
      standardization_params = state_standardization_params
  ) 

  train_df <- scale_res$df
  final_params <- scale_res$standardization_params
   
  train_df$formatted_date <- train_df$observation_date

  return (list(train_df = train_df, standardization_params = final_params))
}


#' Prepare Test Data
#' 
#' Loads the 2018 dataset and standardizes it using the parameters from Training.
#'
prepare_test_data <- function (
    state_covs, 
    obs_covs,
    cov_tif,
    standardization_params, # Must pass params from training
    placeholder_spec_name = "AMCR"
){
  
  # 1. Load Data
  test_filename <- paste0(placeholder_spec_name, "_zf_filtered_region_2018.csv")
  file_path <- file.path("checklist_data","species", placeholder_spec_name, test_filename)
  
  if (!file.exists(file_path)) {
    stop(paste("Test file not found:", file_path))
  }

  test_df_og <- read.delim(file_path, sep = ",", header = TRUE)
  
  # 2. Filter
  test_df_og <- test_df_og[!is.na(test_df_og$duration_minutes),]
  test_df_og <- test_df_og[
    test_df_og$observation_date >= "2018-05-15" & 
    test_df_og$observation_date <= "2018-07-09",
  ]
  test_df <- test_df_og
  
  # 3. Extract Raw Raster Values
  test_env_df <- extract_state_covs(test_df, cov_tif) 
  
  test_df <- inner_join(test_df, test_env_df, by = "checklist_id")

  # 4. Standardize using Training Params
  scale_res <- standardize_ds(
      test_df, 
      obs_covs, 
      state_covs = state_covs,
      standardization_params = standardization_params
  )

  test_df <- scale_res$df
  
  return (test_df) 
}


#' Generate Spatial Structures for Test Data
#' 
#' Creates Voronoi geometries, calculates site areas, and builds the sparse W matrix.
#'
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
  
  # Note: voronoi_clipped_buffers is sourced from model_helpers.R
  test_geoms <- voronoi_clipped_buffers(test_sf_albers, buffer_dist = buffer_m)
  
  # 4. Calculate Site Area (Static)
  # Sum of raster cell areas within the buffer
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
  
  # --- FILTER NA CELLS (New Safety Check) ---
  # Remove cells that fall off the map (NA values) to prevent model errors
  val_col_name <- names(cov_raster_albers)[1]
  test_overlap_df <- test_overlap_df[!is.na(test_overlap_df[[val_col_name]]), ]

  # Calculate Area Weights
  r_res <- terra::res(cov_raster_albers)
  cell_area_km2 <- (r_res[1] / 1000) * (r_res[2] / 1000)
  
  test_overlap_df$w_area <- test_overlap_df$fraction * cell_area_km2

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
  
  return(list(
    test_df = test_df,
    w_matrix = w_matrix
  ))
}