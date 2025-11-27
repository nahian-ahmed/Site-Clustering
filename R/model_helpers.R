library(sf)
library(dplyr)
library(terra)
library(Matrix)
library(tidyr)
library(unmarked)


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


#' Prepare Data for occuN Model
#'
#' Handles the joining of checklists to sites, re-indexing site IDs,
#' pivoting data to wide format, and creating the unmarkedFrame.
prepare_occuN_data <- function(train_data, clustering_df, w_matrix, obs_cov_names, cell_covs) {
  
  # 1. Join checklists to the specific clustering method's site IDs
  comparison_site_lookup <- clustering_df %>%
    dplyr::select(checklist_id, comparison_site = site)
  
  train_data_prepped <- train_data %>%
    dplyr::left_join(comparison_site_lookup, by = "checklist_id") %>%
    dplyr::filter(!is.na(comparison_site)) %>%
    dplyr::select(-site) %>%          
    dplyr::rename(site = comparison_site) 
  
  # 2. Create numeric site index matching the W matrix rows
  M <- nrow(w_matrix)
  site_id_lookup <- data.frame(
    site_char = rownames(w_matrix),  
    site_numeric = 1:M
  )
  
  # Ensure character matching
  train_data_prepped$site <- as.character(train_data_prepped$site)
  site_id_lookup$site_char <- as.character(site_id_lookup$site_char)
  
  # 3. Re-index and create visit IDs
  train_data_prepped <- train_data_prepped %>%
    dplyr::inner_join(site_id_lookup, by = c("site" = "site_char")) %>%
    dplyr::select(-site) %>%
    dplyr::rename(site = site_numeric) %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(visit_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # 4. Pivot Observations (y)
  # Join with 1:M to ensure all sites in W exist in Y, even if no detections
  y_wide <- train_data_prepped %>% 
    tidyr::pivot_wider(
      id_cols = site,
      names_from = visit_id,
      values_from = species_observed,
      values_fill = NA 
    ) %>%
    dplyr::right_join(data.frame(site = 1:M), by = "site") %>%
    dplyr::arrange(site) %>%
    dplyr::select(-site) %>%
    as.matrix()
  
  # 5. Pivot Observation Covariates
  obs_covs_wide <- list()
  for (cov_name in obs_cov_names) {
    obs_covs_wide[[cov_name]] <- train_data_prepped %>% 
      tidyr::pivot_wider(
        id_cols = site,
        names_from = visit_id,
        values_from = dplyr::all_of(cov_name),
        values_fill = NA
      ) %>%
      dplyr::right_join(data.frame(site = 1:M), by = "site") %>%
      dplyr::arrange(site) %>%
      dplyr::select(-site) %>%
      as.matrix()
  }
  
  # 6. Create Unmarked Frame
  umf <- unmarked::unmarkedFrameOccuN(
    y = y_wide,
    obsCovs = obs_covs_wide,
    cellCovs = cell_covs,
    w = w_matrix
  )
  
  return(umf)
}

# #' Fit occuN Model with Random Starts
# #'
# #' repeatedly fits the model to find the global minimum NLL.
# fit_occuN_model <- function(umf, state_formula, obs_formula, n_reps = 30, optimizer = "nlminb") {
  
#   # Combine formulas
#   occuN_formula <- as.formula(paste(
#     paste(deparse(obs_formula), collapse = ""), 
#     paste(deparse(state_formula), collapse = "")
#   ))
  
#   # Determine number of parameters
#   # We fit once briefly or inspect design to get param count, 
#   # or calculate based on covariates. 
#   # Assuming basic formula structure:
#   n_obs_pars <- length(all.vars(obs_formula)) + 1 # +1 intercept
#   n_state_pars <- length(all.vars(state_formula)) + 1
#   n_params <- n_obs_pars + n_state_pars
  
#   best_fm <- NULL
#   min_nll <- Inf
#   fit_successful <- FALSE
  
#   # Loop for random starts
#   for (rep in 1:n_reps) {
#     rand_starts <- runif(n_params, -2, 2) # Slightly tighter bounds usually safer
    
#     fm_rep <- try(unmarked::occuN(
#       formula = occuN_formula,
#       data = umf,
#       starts = rand_starts,
#       se = FALSE, # Optimization speedup
#       method = optimizer
#     ), silent = TRUE)
    
#     if (!inherits(fm_rep, "try-error")) {
#       # Basic check for valid convergence (code 0 or 1 usually ok in R optim)
#       if (fm_rep@negLogLike < min_nll && !is.nan(fm_rep@negLogLike)) {
#         min_nll <- fm_rep@negLogLike
#         best_fm <- fm_rep
#         fit_successful <- TRUE
#       }
#     }
#   }
  
#   if (!fit_successful) return(NULL)
  
#   # Refit best model with SE=TRUE to get Hessians/SEs if needed later, 
#   # or just return the best object found.
#   return(best_fm)
# }



#' Fit occuN Model with Random Starts and Early Stopping
#'
#' Repeatedly fits the model to find the global minimum NLL.
#' Stops early if the same minimum is found multiple times.
fit_occuN_model <- function(umf, state_formula, obs_formula, n_reps = 30, stable_reps = 5, optimizer = "nlminb") {
  
  # Combine formulas
  occuN_formula <- as.formula(paste(
    paste(deparse(obs_formula), collapse = ""), 
    paste(deparse(state_formula), collapse = "")
  ))
  
  # Determine number of parameters to generate valid random starts
  # Assuming basic formula structure (1 intercept + N covariates)
  n_obs_pars <- length(all.vars(obs_formula)) + 1 
  n_state_pars <- length(all.vars(state_formula)) + 1
  n_params <- n_obs_pars + n_state_pars
  
  best_fm <- NULL
  min_nll <- Inf
  fit_successful <- FALSE
  
  # Early stopping counters
  stable_count <- 0
  tolerance <- 0.01  # NLL difference threshold to consider "same result"
  
  # Loop for random starts
  for (rep in 1:n_reps) {
    rand_starts <- runif(n_params, -2, 2) # Random starts [-2, 2]
    
    fm_rep <- try(unmarked::occuN(
      formula = occuN_formula,
      data = umf,
      starts = rand_starts,
      se = FALSE, # Disable SE calculation during search for speed
      method = optimizer
    ), silent = TRUE)
    
    if (!inherits(fm_rep, "try-error")) {
      current_nll <- fm_rep@negLogLike
      
      # Basic check for valid convergence
      if (!is.nan(current_nll)) {
        
        # Case 1: We found a strictly better model
        if (current_nll < min_nll) {
          
          # Check if improvement is significant or just numerical noise
          if (abs(min_nll - current_nll) < tolerance) {
             stable_count <- stable_count + 1
          } else {
             # Significant improvement found: reset stability counter
             stable_count <- 0
          }
          
          min_nll <- current_nll
          best_fm <- fm_rep
          fit_successful <- TRUE
          
        } 
        # Case 2: We hit the existing best minimum again (within tolerance)
        else if (abs(current_nll - min_nll) < tolerance) {
          stable_count <- stable_count + 1
        }
        
        # Early Stopping: If we've hit the same best minimum 3 times, stop.
        if (stable_count >= stable_reps) {
          # cat(sprintf("    Converged early at rep %d\n", rep)) # Optional debug
          break 
        }
      }
    }
  }
  
  if (!fit_successful) return(NULL)
  
  return(best_fm)
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



