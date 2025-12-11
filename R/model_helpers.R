library(sf)
library(dplyr)
library(terra)
library(Matrix)
library(tidyr)
library(unmarked)


########
calculate_weighted_sum <- function(pars, covs_df, intercept = TRUE){

  # --- Vectorized Logic ---
  par_names <- names(pars)
  par_values <- unlist(pars)
  
  if (intercept) {
    intercept_val <- par_values[1]
    cov_names <- par_names[-1]
    cov_values <- par_values[-1]
    
    # Ensure columns are in the exact same order as cov_values
    X_matrix <- as.matrix(covs_df[, cov_names, drop = FALSE])
    
    # Calculate (X %*% B) + intercept
    weighted_sum <- as.numeric(intercept_val + (X_matrix %*% cov_values))
    
  } else {
    cov_names <- par_names
    cov_values <- par_values
    
    X_matrix <- as.matrix(covs_df[, cov_names, drop = FALSE])
    
    weighted_sum <- as.numeric(X_matrix %*% cov_values)
  }
  
  return(weighted_sum)
}

#' Generate Jigsaw-style Site Geometries (Smoothed & Robust)
voronoi_clipped_buffers <- function(points_sf, buffer_dist) {
  
  # 1. Prepare Envelope
  bbox_polygon <- sf::st_as_sfc(sf::st_bbox(points_sf) + buffer_dist * 3)
  
  # 2. Generate Voronoi Tiles for ALL points (Base Partition)
  voronoi_tiles <- sf::st_voronoi(sf::st_union(points_sf), envelope = bbox_polygon, dTolerance = 0.1) %>%
    sf::st_collection_extract(type = "POLYGON") %>%
    sf::st_sf()
  
  # 3. Map Tiles back to Site IDs
  voronoi_w_id <- sf::st_join(voronoi_tiles, points_sf, join = sf::st_contains)
  
  # 4. Define "Exclusive Domain" (Voronoi Territory)
  site_territories <- voronoi_w_id %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    sf::st_make_valid()
  
  # 5. Define "Ideal Domain" (Buffered Convex Hull)
  #    Connects points in a site and adds body
  site_buffers <- points_sf %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(geometry = sf::st_convex_hull(sf::st_union(geometry)), .groups = "drop") %>%
    sf::st_buffer(dist = buffer_dist) %>%
    sf::st_make_valid()
  
  # 6. Intersect to get Final Shapes
  site_territories <- dplyr::rename(site_territories, site_t = site)
  site_buffers <- dplyr::rename(site_buffers, site_b = site)
  

  # Define attributes as constant before intersection
  sf::st_agr(site_territories) <- "constant"
  sf::st_agr(site_buffers) <- "constant"

  # Calculate intersection (The Cookie Cutter)
  intersections <- sf::st_intersection(site_territories, site_buffers)
  
  final_geoms <- intersections %>%
    dplyr::filter(site_t == site_b) %>%
    dplyr::rename(site = site_t) %>%
    dplyr::select(site, geometry) %>%
    sf::st_make_valid() 
  
  # 7. [NEW] Smoothing Step
  #    removes sharp Voronoi corners while keeping general shape
  #    dTolerance = 10m provides a subtle smoothing effect
  final_geoms <- sf::st_simplify(final_geoms, dTolerance = 10, preserveTopology = TRUE)

  # 8. [FIXED] Robust Collection Extraction
  #    Only run extraction if the geometry is actually a collection.
  #    This prevents the "x is already of type POLYGON" warning.
  geom_type <- sf::st_geometry_type(final_geoms, by_geometry = FALSE)
  if (inherits(geom_type, "GEOMETRYCOLLECTION") || any(geom_type == "GEOMETRYCOLLECTION")) {
      final_geoms <- sf::st_collection_extract(final_geoms, "POLYGON")
  }

  return(final_geoms)
}

create_site_geometries <- function(site_data_sf, reference_raster, buffer_m = 15, method_name = NULL, area_unit = "m") {
  
  # --- 1. Get CRS from the Reference Raster ---
  # This ensures we are working in the exact same projection (Albers 100m)
  # as the simulation grid.
  target_crs <- sf::st_crs(terra::crs(reference_raster))
  
  # Project Points to Target CRS (Albers)
  points_sf <- sf::st_as_sf(
    site_data_sf,
    coords = c("longitude", "latitude"),
    crs = 4326, # WGS84
    remove = FALSE 
  )
  
  points_proj <- sf::st_transform(points_sf, crs = target_crs)
  
  # --- 2. Create Site Geometries ---
  is_kmsq <- !is.null(method_name) && grepl("kmSq", method_name, ignore.case = TRUE)
  
  if (is_kmsq) {
    # --- A. SQUARE GRID GEOMETRIES (kmSq) ---
    parts <- unlist(strsplit(method_name, "-"))
    rad_m <- NA
    
    if (parts[1] == "kmSq") {
      rad_m <- as.numeric(parts[2])
    } else {
      area_km <- as.numeric(parts[1])
      rad_m <- as.integer(sqrt(area_km) * 1000)
    }
    
    bbox <- sf::st_bbox(points_proj)
    
    # Create grid covering the extent
    x_seq <- seq(from = bbox["xmin"] - (rad_m * 5), to = bbox["xmax"] + (rad_m * 5), by = rad_m)
    y_seq <- seq(from = bbox["ymin"] - (rad_m * 5), to = bbox["ymax"] + (rad_m * 5), by = rad_m)
    
    grid_points <- expand.grid(x = x_seq, y = y_seq)
    
    # Use target_crs here
    grid_points_sf <- sf::st_as_sf(grid_points, coords = c("x", "y"), crs = target_crs)
    
    grid_geom <- sf::st_make_grid(grid_points_sf, cellsize = rad_m, square = TRUE)
    grid_sf <- sf::st_sf(geometry = grid_geom)
    grid_sf$site <- as.character(seq_len(nrow(grid_sf)))
    
    present_sites <- unique(as.character(site_data_sf$site))
    site_geoms_sf <- grid_sf[grid_sf$site %in% present_sites, ]
    
  } else {
    # --- B. POINT-BASED VORONOI CLIP ---
    # Pass the projected points (meters) to the Voronoi function
    site_geoms_sf <- voronoi_clipped_buffers(points_proj, buffer_m)
  }
  
  # --- 3. Create the 'w' (overlap weight) matrix ---
  site_geoms_sf$site <- as.character(site_geoms_sf$site)
  
  site_vect <- terra::vect(site_geoms_sf)
  # No need to project site_vect again if we built it using target_crs, 
  # but keeping it is a safe double-check.
  site_vect_proj <- terra::project(site_vect, terra::crs(reference_raster))
  
  overlap_df <- terra::extract(
    reference_raster[[1]], 
    site_vect_proj, 
    cells = TRUE, 
    exact = TRUE, 
    ID = TRUE
  )
  
  overlap_df$w_area <- overlap_df$fraction
  
  n_sites <- nrow(site_geoms_sf)
  n_cells <- terra::ncell(reference_raster)
  
  w <- Matrix::sparseMatrix(
    i = overlap_df$ID,
    j = overlap_df$cell,
    x = overlap_df$w_area,
    dims = c(n_sites, n_cells),
    dimnames = list(site_geoms_sf$site, NULL) 
  )
  
  attr(site_geoms_sf, "w_matrix") <- w
  
  return(site_geoms_sf)
}

split_disjoint_sites <- function(site_geoms_sf, method_name, min_area_threshold = 100) {
  
  # 1. SKIP "SAFE" METHODS
  # kmSq (grids), rounded (coords), lat-long/UL (points) are always contiguous.
  if (grepl("kmSq", method_name, ignore.case = TRUE) ||
      grepl("rounded", method_name, ignore.case = TRUE) ||
      grepl("lat-long", method_name, ignore.case = TRUE) ||
      grepl("UL", method_name, ignore.case = TRUE) ||
      grepl("SVS", method_name, ignore.case = TRUE)) {
    return(site_geoms_sf)
  }
  
  cat(sprintf("    [Post-Process] Checking %s for disjoint geometries...\n", method_name))
  
  # 2. CAST MULTIPOLYGONS TO POLYGONS
  # This explodes disjoint parts into separate rows
  # Suppress warnings about attributes being constant
  sites_split <- suppressWarnings(sf::st_cast(site_geoms_sf, "POLYGON", warn = FALSE))
  
  # If no splitting happened, return original (save computation)
  if (nrow(sites_split) == nrow(site_geoms_sf)) {
    return(site_geoms_sf)
  }
  
  # 3. FILTER SLIVERS
  # Remove tiny artifacts created by the Voronoi clipping process
  sites_split$area_tmp <- as.numeric(sf::st_area(sites_split))
  sites_clean <- sites_split[sites_split$area_tmp > min_area_threshold, ]
  sites_clean$area_tmp <- NULL # Cleanup
  
  # 4. RE-GENERATE SITE IDs
  # Append a suffix to the original site ID to create unique IDs for the new parts
  # e.g., Site "10" splits into "10_1" and "10_2"
  # We group by the OLD site ID to generate the sequence
  sites_clean <- sites_clean %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(
      part_id = dplyr::row_number(),
      new_site_id = paste0(site, "_", part_id)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-site, -part_id) %>%
    dplyr::rename(site = new_site_id)
  
  cat(sprintf("    [Post-Process] Split disjoint sites. Count changed from %d to %d.\n", 
              nrow(site_geoms_sf), nrow(sites_clean)))
  
  return(sites_clean)
}

#' Prepare Data for occuN Model
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


#' Fit occuN Model with Random Starts and Early Stopping
fit_occuN_model <- function(umf, state_formula, obs_formula, n_reps = 30, 
                            stable_reps = 10, optimizer = "nlminb", 
                            lower = -Inf, upper = Inf) { # <--- 1. Add arguments here
  
  occuN_formula <- as.formula(paste(
    paste(deparse(obs_formula), collapse = ""), 
    paste(deparse(state_formula), collapse = "")
  ))
  
  # Note: This parameter counting method assumes simple formulas (intercept + linear terms).
  # If you use interactions (*) or factors, model.matrix would be safer to count columns.
  n_obs_pars <- length(all.vars(obs_formula)) + 1 
  n_state_pars <- length(all.vars(state_formula)) + 1
  n_params <- n_obs_pars + n_state_pars
  
  best_fm <- NULL
  min_nll <- Inf
  fit_successful <- FALSE
  
  stable_count <- 0
  tolerance <- 0.01 
  
  for (rep in 1:n_reps) {
    # Generate random starts
    # Note: If your bounds are tight (e.g., lower=0), 
    # make sure your starts respect them!
    rand_starts <- runif(n_params, -2, 2)
    
    # If bounds are provided, clamp the starts to be within bounds
    # to avoid immediate rejection by the optimizer
    if(any(is.finite(lower)) || any(is.finite(upper))) {
       # Simple clamping logic
       safe_lower <- ifelse(is.finite(lower), lower + 0.01, -2)
       safe_upper <- ifelse(is.finite(upper), upper - 0.01, 2)
       # Ensure we don't cross (if lower > -2, use lower, etc)
       # This is a basic safety check, adjust logic if specific bounds are needed
       rand_starts <- pmax(pmin(rand_starts, safe_upper), safe_lower)
    }

    fm_rep <- try(unmarked::occuN(
      formula = occuN_formula,
      data = umf,
      starts = rand_starts,
      se = FALSE, 
      method = optimizer,
      lower = lower,
      upper = upper
    ), silent = TRUE)
    
    if (!inherits(fm_rep, "try-error")) {
      current_nll <- fm_rep@negLogLike
      
      if (!is.nan(current_nll)) {
        if (current_nll < min_nll) {
          if (abs(min_nll - current_nll) < tolerance) {
             stable_count <- stable_count + 1
          } else {
             stable_count <- 0
          }
          min_nll <- current_nll
          best_fm <- fm_rep
          fit_successful <- TRUE
        } 
        else if (abs(current_nll - min_nll) < tolerance) {
          stable_count <- stable_count + 1
        }
        
        if (stable_count >= stable_reps) {
          break 
        }
      }
    }
  }
  
  if (!fit_successful) return(NULL)
  return(best_fm)
}


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

  region$occ_prob <- rje::expit(weighted_sum)

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