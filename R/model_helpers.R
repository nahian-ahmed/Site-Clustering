library(sf)
library(dplyr)
library(terra)
library(Matrix)

##########
# site closure for Occ Model:
#   1. constant site covariates
#   2. no false positives (detected only if occupied)
##########
enforceClosure <- function(sites_df, occ_cov_list, sites_list){
  j<-1
  closed_df <- NA

  for(eBird_site in sites_list){
    
  
    checklists_at_site <- sites_df[sites_df$site == eBird_site,]
    
    for(occCov_i in occ_cov_list){
      checklists_at_site[occCov_i] <- mean(checklists_at_site[[occCov_i]])
    
    }

    
    if(j==1){
      closed_df = checklists_at_site
    } else {
      closed_df = rbind(closed_df, checklists_at_site)
    }
    j = j+1
  
  }
  return(closed_df)
}


#######
# calculates the occupancy model from a given dataset
# containing checklists and sites
#######
# 1. enforces closure
# 2. formats it w/r/t eBird data
# 3. runs through occupancy model
#######
calcOccModel <- function(df, occ_covs, det_covs, skip_closure=FALSE){

  sites_occ <- subset(df, !duplicated(site))$site

  
  closed_df <- df
  if(!skip_closure){
    closed_df <- enforceClosure(df, occ_covs, sites_occ)
  } 
  
  
  umf_AUK <- auk::format_unmarked_occu(
    closed_df,
    site_id = "site",
    response = "species_observed",
    site_covs = occ_covs,
    obs_covs = det_covs
  )
  
  det_cov_str <- paste("", paste(det_covs, collapse="+"), sep=" ~ ")
  occ_cov_str <- paste("", paste(occ_covs, collapse="+"), sep=" ~ ")
  
  species_formula <- paste(det_cov_str, occ_cov_str, sep = " ")
  species_formula <- as.formula(species_formula)
  
  occ_um <- unmarked::formatWide(umf_AUK, type = "unmarkedFrameOccu")

  og_syn_gen_form <- unmarked::occu(formula = species_formula, data = occ_um)
  
  return(og_syn_gen_form)

}
#######





########
predict_sdm_map <- function(occ_pars, region, intercept = TRUE){
  
  valid_boundary <- terra::vect("occupancy_feature_raster/boundary/boundary.shp")
  crs(valid_boundary) <- crs(region)
  region <- terra::crop(region, valid_boundary, mask = TRUE)
  
  weighted_sum <- 0
  par_idx <- 1
  n_pars <- length(occ_pars)

  if(intercept){
    weighted_sum <- weighted_sum + occ_pars[[1]]
    par_idx <- par_idx + 1
  }
  while(par_idx <= n_pars){
    weighted_sum <- weighted_sum + (occ_pars[[par_idx]] * subset(region, names(occ_pars)[[par_idx]]))
    par_idx <- par_idx + 1
  }

  region$occ_prob <- expit(weighted_sum)


  return(region)
}
########

########
get_occu_map_diff <- function(occu_map_gt, occu_map) {

  occu_m_gt <- occu_map_gt[["occu"]]
  occu_m <- occu_map[["occu"]]

  difference <- occu_m - occu_m_gt
  percentage_difference <- (abs(occu_m - occu_m_gt) / occu_m_gt) * 100

  # Create a raster stack with two bands
  occu_map_diff <- rast(nrows = nrow(occu_m_gt), ncols = ncol(occu_m_gt), 
              ext = ext(occu_m_gt), crs = crs(occu_m_gt), nlyrs = 2)
  
  values(occu_map_diff)[,1] <- values(difference)
  values(occu_map_diff)[,2] <- values(percentage_difference)

  names(occu_map_diff) <- c("difference", "percentage_difference")

  return(occu_map_diff)
}
########


########
calculate_weighted_sum <- function(pars, covs_df, intercept = TRUE){

  # --- NEW Vectorized Logic ---
  # This is much faster than the original R loop.
  
  # 1. Get parameter names and values
  # The order is guaranteed to be consistent
  par_names <- names(pars)
  par_values <- unlist(pars)
  
  if (intercept) {
    # 2a. Separate intercept from covariates
    intercept_val <- par_values[1]
    cov_names <- par_names[-1]
    cov_values <- par_values[-1]
    
    # 3a. Get the covariate matrix (X)
    # Ensure columns are in the *exact* same order as cov_values
    X_matrix <- as.matrix(covs_df[, cov_names, drop = FALSE])
    
    # 4a. Calculate (X %*% B) + intercept
    # as.numeric() ensures it's a simple vector, not a 1-column matrix
    weighted_sum <- as.numeric(intercept_val + (X_matrix %*% cov_values))
    
  } else {
    # 2b. No intercept, use all parameters
    cov_names <- par_names
    cov_values <- par_values
    
    # 3b. Get the covariate matrix (X)
    X_matrix <- as.matrix(covs_df[, cov_names, drop = FALSE])
    
    # 4b. Calculate (X %*% B)
    weighted_sum <- as.numeric(X_matrix %*% cov_values)
  }
  
  return(weighted_sum)
}


########
get_parameters <- function(df, i, occ_covs, det_covs, occ_intercept = TRUE, det_intercept = TRUE){

  occ_par_list <- list()
  if (occ_intercept){
    occ_par_list[["occ_intercept"]] <- df[i, "occ_intercept"] 
  }
  for (occ_cov in occ_covs){
    occ_par_list[[occ_cov]] <- df[i, occ_cov]

  }

  det_par_list <- list()
  if (det_intercept){
    det_par_list[["det_intercept"]] <- df[i, "det_intercept"] 
  }
  for (det_cov in det_covs){
    det_par_list[[det_cov]] <- df[i, det_cov]

  }

 
  return (list(occ_par_list = occ_par_list, det_par_list = det_par_list))
}
########


create_site_geometries <- function(site_data_sf, reference_raster, buffer_m = 15, method_name = NULL, area_unit = "m") {
  
  # --- 1. Project Points to Albers (meters) for Geometry Creation ---
  albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  
  points_sf <- sf::st_as_sf(
    site_data_sf,
    coords = c("longitude", "latitude"),
    crs = "+proj=longlat +datum=WGS84",
    remove = FALSE 
  )
  
  points_albers <- sf::st_transform(points_sf, crs = albers_crs_str)
  
  
  # --- 2. Create Site Geometries (in Albers) ---
  is_kmsq <- !is.null(method_name) && grepl("kmSq", method_name, ignore.case = TRUE)
  
  if (is_kmsq) {
    # --- A. SQUARE GRID GEOMETRIES (kmSq) ---
    parts <- unlist(strsplit(method_name, "-"))
    rad_m <- NA
    
    if (parts[1] == "kmSq") {
      # Format: "kmSq-1000"
      rad_m <- as.numeric(parts[2])
    } else {
      # Format: "1-kmSq" or "0.125-kmSq"
      area_km <- as.numeric(parts[1])
      
      # FIX: Use as.integer to match R/clustering_helpers.R logic
      # This ensures the grid used for plotting aligns exactly with the grid used for clustering
      rad_m <- as.integer(sqrt(area_km) * 1000)
    }
    
    # 1. Get bounding box 
    bbox <- sf::st_bbox(points_albers)
    
    # 2. Define grid extent (using the exact same expansion logic as kmsq.R)
    # Note: If site_data_sf is a subset, this bbox might differ from the original clustering bbox.
    # Ideally, you should pass the full dataset's bbox, but matching rad_m fixes the main drift.
    x_seq <- seq(from = bbox["xmin"] - (rad_m * 5), to = bbox["xmax"] + (rad_m * 5), by = rad_m)
    y_seq <- seq(from = bbox["ymin"] - (rad_m * 5), to = bbox["ymax"] + (rad_m * 5), by = rad_m)
    
    # 3. Create grid points
    grid_points <- expand.grid(x = x_seq, y = y_seq)
    grid_points_sf <- sf::st_as_sf(grid_points, coords = c("x", "y"), crs = albers_crs_str)
    
    # 4. Create square polygons
    grid_geom <- sf::st_make_grid(grid_points_sf, cellsize = rad_m, square = TRUE)
    grid_sf <- sf::st_sf(geometry = grid_geom)
    
    # 5. Assign IDs 
    grid_sf$site <- as.character(seq_len(nrow(grid_sf)))
    
    # 6. Filter to keep only sites present in the data
    present_sites <- unique(as.character(site_data_sf$site))
    site_geoms_sf <- grid_sf[grid_sf$site %in% present_sites, ]
    
  } else {
    # --- B. CONVEX HULL + BUFFER ---
    site_geoms_sf <- points_albers %>%
      group_by(site) %>%
      summarise(
        geometry = sf::st_buffer(
            sf::st_convex_hull(sf::st_union(geometry)),
            dist = buffer_m
        ),
        .groups = "drop" 
      )
  }
  
  # --- 3. Create the 'w' (overlap weight) matrix ---
  site_vect <- terra::vect(site_geoms_sf)
  site_vect_proj <- terra::project(site_vect, terra::crs(reference_raster))
  
  overlap_df <- terra::extract(
    reference_raster[[1]], 
    site_vect_proj, 
    cells = TRUE, 
    exact = TRUE, 
    ID = TRUE
  )
  
  cell_areas <- terra::cellSize(reference_raster, unit = "m")
  overlap_cell_areas <- terra::extract(cell_areas, overlap_df$cell)
  
  overlap_df$w_area <- overlap_df$fraction * overlap_cell_areas[,1]
  
  if (area_unit == "km") {
    overlap_df$w_area <- overlap_df$w_area / 1e6 
  }
  
  n_sites <- nrow(site_geoms_sf)
  n_cells <- terra::ncell(reference_raster)
  
  w <- Matrix::sparseMatrix(
    i = overlap_df$ID,
    j = overlap_df$cell,
    x = overlap_df$w_area,
    dims = c(n_sites, n_cells),
    dimnames = list(as.character(site_geoms_sf$site), NULL) 
  )
  
  attr(site_geoms_sf, "w_matrix") <- w
  
  return(site_geoms_sf)
}