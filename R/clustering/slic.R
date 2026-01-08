###################################################################
# SLIC Clustering Helper (slic.R)
#
# Implements Simple Linear Iterative Clustering (SLIC) for 
# spatially-contiguous site generation, adapted for eBird data.
#
# Generated: January 08, 2026
###################################################################

library(terra)
library(sf)
library(dplyr)
library(stats)

# -----------------------------------------------------------------
# 1. Generate Initial Seeds
# -----------------------------------------------------------------

#' Generate Initial Seeds for SLIC
#' 
#' Creates a regular grid of seed points based on zeta, then filters 
#' them to keep only seeds within a specified buffer of existing checklists.
#' 
#' @param r Template raster (for extent and CRS)
#' @param zeta Target number of seeds (determines grid spacing)
#' @param checklists_sf sf object containing checklist locations
#' @param buffer_dist_m Distance buffer to valid seeds (default 50km)
#' 
#' @return sf object of valid seed points
get_slic_seeds <- function(r, zeta, checklists_sf, buffer_dist_m = 50000) {
  
  # --- A. Define Regular Grid ---
  ext_r <- ext(r)
  width <- ext_r$xmax - ext_r$xmin
  height <- ext_r$ymax - ext_r$ymin
  area_m2 <- width * height
  
  # Calculate grid interval S = sqrt(Area / zeta)
  S <- sqrt(area_m2 / zeta)
  
  # Generate grid coordinates centered in cells
  x_coords <- seq(ext_r$xmin + S/2, ext_r$xmax, by = S)
  y_coords <- seq(ext_r$ymin + S/2, ext_r$ymax, by = S)
  
  grid_df <- expand.grid(x = x_coords, y = y_coords)
  seeds_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = crs(r))
  
  # --- B. Filter Seeds Near Data ---
  # Ensure checklists are in the same CRS
  checklists_sf <- st_transform(checklists_sf, crs(r))
  
  # Use efficient spatial indexing to find seeds within buffer
  # sparse = TRUE returns indices; if length > 0, the seed is valid
  is_near <- st_is_within_distance(seeds_sf, checklists_sf, dist = buffer_dist_m, sparse = TRUE)
  keep_indices <- lengths(is_near) > 0
  
  valid_seeds <- seeds_sf[keep_indices, ]
  
  return(valid_seeds)
}


# -----------------------------------------------------------------
# 2. Perform SLIC Algorithm on Raster Pixels
# -----------------------------------------------------------------

#' Run SLIC Algorithm
#' 
#' Assigns every valid raster pixel to a cluster using a spatially-weighted
#' K-means approach (SLIC adaptation).
#' 
#' @param cov_raster Raster stack of covariates (standardized)
#' @param seeds_sf sf object of initial seed locations
#' @param eta Compactness parameter (controls weight of spatial vs feature distance)
#' @param zeta Initial number of seeds (used to normalize spatial distance)
#' 
#' @return Dataframe of valid raster cells with columns: x, y, [covs], site
perform_slic_clustering <- function(cov_raster, seeds_sf, eta, zeta) {
  
  # 1. Convert Raster to Dataframe (Valid Pixels Only)
  df_pixels <- as.data.frame(cov_raster, xy = TRUE, na.rm = TRUE)
  
  if (nrow(df_pixels) == 0) stop("Raster contains no valid pixels.")
  
  # Identify covariate names from the dataframe itself (to handle sanitization)
  state_cov_names <- setdiff(names(df_pixels), c("x", "y"))
  
  # 2. Extract Features for Initial Seeds
  seed_vals <- terra::extract(cov_raster, seeds_sf)
  
  # --- ROBUST NAME ALIGNMENT START ---
  # terra::extract() often returns an 'ID' column and might not sanitize names 
  # exactly the same way as.data.frame() does. We fix this here.
  
  # Remove 'ID' column if present
  if("ID" %in% names(seed_vals)) {
    seed_vals <- seed_vals[, setdiff(names(seed_vals), "ID"), drop = FALSE]
  }
  
  # Ensure column count matches
  if(ncol(seed_vals) != length(state_cov_names)) {
    stop(paste("Mismatch in column counts: Raster DF has", length(state_cov_names), 
               "covariates, but extracted seeds have", ncol(seed_vals)))
  }
  
  # Force names to match df_pixels (Order is preserved by terra)
  names(seed_vals) <- state_cov_names
  # --- ROBUST NAME ALIGNMENT END ---
  
  seed_coords <- st_coordinates(seeds_sf)
  
  # Combine features and coords
  seed_data <- cbind(seed_vals, as.data.frame(seed_coords))
  
  # Remove seeds that fall on NA pixels (outside mask)
  seed_data <- na.omit(seed_data)
  
  if (nrow(seed_data) == 0) stop("No valid seeds found on the raster mask.")
  
  # 3. Scaling for SLIC Distance Metric
  # Distance D^2 = d_feat^2 + (eta/S)^2 * d_xy^2
  
  ex <- ext(cov_raster)
  area_domain <- (ex$xmax - ex$xmin) * (ex$ymax - ex$ymin)
  S <- sqrt(area_domain / zeta)
  
  spatial_scale <- eta / S
  
  # 4. Prepare Matrices for Lloyd's Algorithm (K-Means)
  
  # -- Data Matrix (Pixels) --
  feat_cols <- as.matrix(df_pixels[, state_cov_names, drop = FALSE])
  xy_cols <- as.matrix(df_pixels[, c("x", "y")]) * spatial_scale
  data_matrix <- cbind(feat_cols, xy_cols)
  
  # -- Centers Matrix (Seeds) --
  # Must match column order of data_matrix
  seed_feat <- as.matrix(seed_data[, state_cov_names, drop = FALSE])
  seed_xy <- as.matrix(seed_data[, c("x", "y")]) * spatial_scale
  centers_matrix <- cbind(seed_feat, seed_xy)
  
  # 5. Run Clustering
  # Using iter.max=20 is standard for SLIC as it converges quickly.
  # algorithm="Lloyd" ensures standard iterative refinement.
  km_res <- kmeans(data_matrix, centers = centers_matrix, iter.max = 20, algorithm = "Lloyd")
  
  # 6. Attach Labels
  df_pixels$site <- km_res$cluster
  
  return(df_pixels)
}


# -----------------------------------------------------------------
# 3. Main Wrapper Function
# -----------------------------------------------------------------

#' SLIC Sites
#' 
#' Main entry point for the pipeline. Generates sites and assigns checklists.
#' 
#' @param checklists Dataframe of observations (must have latitude, longitude)
#' @param state_covs Vector of covariate names
#' @param cov_raster terra SpatRaster object (must be standardized)
#' @param eta Compactness
#' @param zeta Number of initial seeds
#' 
#' @return checklists dataframe with a new 'site' column
slicSites <- function(checklists, state_covs, cov_raster, eta, zeta) {
  
  # Convert checklists to SF for spatial operations
  checklists_sf <- st_as_sf(checklists, coords = c("longitude", "latitude"), crs = 4326)
  
  # Transform checklists to match Raster CRS (usually Albers)
  checklists_sf <- st_transform(checklists_sf, crs(cov_raster))
  
  # 1. Generate Filtered Seeds
  # buffer_dist_m can be tuned; 50km covers most eBird data gaps efficiently
  seeds <- get_slic_seeds(cov_raster, zeta, checklists_sf, buffer_dist_m = 50000)
  
  if (nrow(seeds) == 0) {
    warning("No seeds found near checklists. Assigning all to site -1.")
    checklists$site <- -1
    return(checklists)
  }
  
  # 2. Run SLIC on the Landscape
  # Returns the raster pixels with site labels
  pixel_labels_df <- perform_slic_clustering(cov_raster, seeds, eta, zeta)
  
  # 3. Map Pixel Labels to Checklists
  # Create a temporary classification raster from the result
  label_raster <- rast(cov_raster, nlyrs = 1)
  names(label_raster) <- "site"
  
  # Fill raster efficiently using cell IDs
  cell_ids <- cellFromXY(label_raster, pixel_labels_df[, c("x", "y")])
  label_raster[cell_ids] <- pixel_labels_df$site
  
  # Extract site ID for each checklist location
  extracted_sites <- terra::extract(label_raster, checklists_sf)
  
  # 4. Attach to original dataframe
  checklists$site <- extracted_sites$site
  
  # Handle checklists falling outside valid raster pixels
  checklists$site[is.na(checklists$site)] <- -1
  
  return(checklists)
}