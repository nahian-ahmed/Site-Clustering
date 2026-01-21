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


#' Create Site Geometries (Updated)
#' 
#' Returns SF object of site boundaries. NO LONGER returns w_matrix.
#'
create_site_geometries <- function(site_data_sf, reference_raster, buffer_m = 15, method_name = NULL) {
  
  target_crs <- sf::st_crs(terra::crs(reference_raster))
  
  # Project Points
  points_sf <- sf::st_as_sf(
    site_data_sf,
    coords = c("longitude", "latitude"),
    crs = 4326, 
    remove = FALSE 
  )
  points_proj <- sf::st_transform(points_sf, crs = target_crs)
  
  # Check Method Type
  is_kmsq <- !is.null(method_name) && grepl("kmSq", method_name, ignore.case = TRUE)
  is_slic <- !is.null(method_name) && grepl("SLIC", method_name, ignore.case = TRUE)
  
  if (is_kmsq) {
    # --- Grid Logic ---
    parts <- unlist(strsplit(method_name, "-"))
    if (parts[1] == "kmSq") {
      rad_m <- as.numeric(parts[2])
    } else {
      area_km <- as.numeric(parts[1])
      rad_m <- as.integer(sqrt(area_km) * 1000)
    }
    
    bbox <- sf::st_bbox(points_proj)
    x_seq <- seq(from = bbox["xmin"] - (rad_m * 5), to = bbox["xmax"] + (rad_m * 5), by = rad_m)
    y_seq <- seq(from = bbox["ymin"] - (rad_m * 5), to = bbox["ymax"] + (rad_m * 5), by = rad_m)
    
    grid_points <- expand.grid(x = x_seq, y = y_seq)
    grid_points_sf <- sf::st_as_sf(grid_points, coords = c("x", "y"), crs = target_crs)
    grid_geom <- sf::st_make_grid(grid_points_sf, cellsize = rad_m, square = TRUE)
    
    grid_sf <- sf::st_sf(geometry = grid_geom)
    grid_sf$site <- as.character(seq_len(nrow(grid_sf)))
    
    present_sites <- unique(as.character(site_data_sf$site))
    site_geoms_sf <- grid_sf[grid_sf$site %in% present_sites, ]
    
  } else if (is_slic) {
    # --- SLIC Logic (Optimized) ---
    parts <- unlist(strsplit(method_name, "-"))
    eta <- as.numeric(parts[2])
    zeta <- as.numeric(parts[3])
    
    # 1. Regenerate Seeds
    seeds <- get_slic_seeds(reference_raster, zeta, points_proj, buffer_dist_m = 50000)
    
    # 2. Run SLIC
    pixel_df <- perform_slic_clustering(reference_raster, seeds, eta, zeta)
    
    # 3. Create Raster
    r_temp <- terra::rast(reference_raster, nlyrs = 1)
    cells <- terra::cellFromXY(r_temp, as.matrix(pixel_df[, c("x", "y")]))
    r_temp[cells] <- pixel_df$site
    
    # --- OPTIMIZATION START ---
    # Identify sites present in the data
    present_sites <- unique(as.character(site_data_sf$site))
    present_sites_num <- as.numeric(present_sites)
    
    # Mask the raster: Set any cell NOT in a present site to NA
    # This prevents creating polygons for empty areas
    r_temp[!r_temp %in% present_sites_num] <- NA
    # --- OPTIMIZATION END ---

    # 4. Convert ONLY valid clusters to Polygons
    polys <- terra::as.polygons(r_temp, dissolve = TRUE)
    site_geoms_sf <- sf::st_as_sf(polys)
    
    # Rename and Final Filter
    colnames(site_geoms_sf)[1] <- "site"
    site_geoms_sf$site <- as.character(site_geoms_sf$site)
    
    # (Double check filter just in case)
    site_geoms_sf <- site_geoms_sf[site_geoms_sf$site %in% present_sites, ]
  } else {
    # --- Voronoi Logic ---
    site_geoms_sf <- voronoi_clipped_buffers(points_proj, buffer_m)
  }
  
  site_geoms_sf$site <- as.character(site_geoms_sf$site)
  
  return(site_geoms_sf)
}


disjoint_site_geometries <- function(site_geoms_sf, point_data_df, crs_points = 4326) {
  
  # 1. Force Cast to MULTIPOLYGON first, then POLYGON
  sites_split <- site_geoms_sf %>%
    sf::st_make_valid() %>%
    sf::st_cast("MULTIPOLYGON") %>% 
    sf::st_cast("POLYGON", warn = FALSE)
  
  # Check if splitting occurred
  if (nrow(sites_split) == nrow(site_geoms_sf)) {
    message("  - No disjoint geometries found. Returning original sites.")
    return(list(geoms = site_geoms_sf, data = point_data_df))
  }
  
  message(sprintf("  - Splitting detected: Increased sites from %d to %d", 
                  nrow(site_geoms_sf), nrow(sites_split)))
  
  # 2. Create Unique IDs for New Parts
  sites_split <- sites_split %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(
      sub_id = dplyr::row_number(),
      new_site_id = paste0(site, "_", sub_id)
    ) %>%
    dplyr::ungroup()
  
  # 3. Robust Point Reassignment
  points_sf <- sf::st_as_sf(
    point_data_df, 
    coords = c("longitude", "latitude"), 
    crs = crs_points, 
    remove = FALSE
  ) %>%
    sf::st_transform(sf::st_crs(sites_split))
  
  # Spatial join to find the NEW polygon for each point
  join_res <- sf::st_join(points_sf, sites_split["new_site_id"], join = sf::st_intersects, left = TRUE)
  
  # Deduplicate if points touch boundaries
  if ("checklist_id" %in% names(join_res)) {
     join_res <- join_res[!duplicated(join_res$checklist_id), ]
  } else {
     join_res <- join_res[!duplicated(geometry), ]
  }
  
  # 4. Update the Dataframe
  point_data_updated <- point_data_df
  
  if ("checklist_id" %in% names(point_data_updated)) {
    match_idx <- match(point_data_updated$checklist_id, join_res$checklist_id)
    new_ids <- join_res$new_site_id[match_idx]
  } else {
    new_ids <- join_res$new_site_id
  }
  
  point_data_updated$site <- ifelse(!is.na(new_ids), new_ids, as.character(point_data_updated$site))
  
  # 5. Finalize Geometries (FILTERING ADDED)
  # Identify which new sites actually contain points
  occupied_new_ids <- unique(na.omit(new_ids))
  
  sites_final <- sites_split %>%
    dplyr::filter(new_site_id %in% occupied_new_ids) %>% # <--- DROPS EMPTY ISLANDS
    dplyr::select(-site, -sub_id) %>%
    dplyr::rename(site = new_site_id) %>%
    dplyr::select(site, geometry)
  
  message(sprintf("  - Filtered empty islands: Reduced to %d occupied sites.", nrow(sites_final)))

  return(list(geoms = sites_final, data = point_data_updated))
}

generate_overlap_matrix <- function(site_geoms_sf, reference_raster) {
  
  # Check if we can use the Fast Raster Path (for SLIC)
  # If the geometry looks like it came from a raster (perfect fit), we could optimize.
  # But the safest way is to check if we can reconstruct it or just optimize the extract.
  
  # Standard Approach (Optimized for NA removal)
  site_geoms_sf$site <- as.character(site_geoms_sf$site)
  site_vect <- terra::vect(site_geoms_sf)
  
  # Check CRS
  if (!identical(crs(site_vect), crs(reference_raster))) {
    site_vect <- terra::project(site_vect, terra::crs(reference_raster))
  }

  # --- OPTIMIZATION FOR SLIC/GRID POLYGONS ---
  # terra::extract on polygons is slow. 
  # If we are memory constrained, we can process in chunks or using 'cells' method carefully.
  # But for SLIC, exact=TRUE is heavy. 
  # Since SLIC aligns with pixels, exact=FALSE (default) is sufficient and much faster
  # because a pixel is either fully in the site or not.
  
  overlap_df <- terra::extract(
    reference_raster[[1]], 
    site_vect, 
    cells = TRUE, 
    exact = FALSE, # Change to FALSE for SLIC/Grid (huge speedup)
    ID = TRUE
  )
  
  # Remove NA cells (water/outside boundary)
  val_col_name <- names(reference_raster)[1]
  overlap_df <- overlap_df[!is.na(overlap_df[[val_col_name]]), ]
  
  # Calculate Area
  r_res <- terra::res(reference_raster)
  cell_area_km2 <- (r_res[1] / 1000) * (r_res[2] / 1000)
  
  # If exact=FALSE, we assume fraction is 1.0 for included cells
  if (!("fraction" %in% names(overlap_df))) {
    overlap_df$fraction <- 1.0
  }
  
  overlap_df$w_area <- overlap_df$fraction * cell_area_km2
  
  # Map IDs back to Site Names
  # extract returns ID as 1,2,3... corresponding to row order of site_vect
  site_indices <- overlap_df$ID
  site_names <- site_geoms_sf$site[site_indices]
  
  # Build Matrix
  n_sites <- nrow(site_geoms_sf)
  n_cells <- terra::ncell(reference_raster)
  
  # We need to map site_names to row indices in the final matrix
  # Let's trust the input order of site_geoms_sf for the rows
  
  w <- Matrix::sparseMatrix(
    i = overlap_df$ID,     # Row index (matches site_geoms_sf order)
    j = overlap_df$cell,   # Column index (raster cell ID)
    x = overlap_df$w_area,
    dims = c(n_sites, n_cells),
    dimnames = list(site_geoms_sf$site, NULL) 
  )
  
  return(w)
}

#' Prepare Data for occuN Model
#' 
#' Updated to automatically filter training data to match the provided clustering_df.
#' This prevents false-alarm warnings when using filtering methods like 1to10.
prepare_occuN_data <- function(train_data, clustering_df, w_matrix, obs_cov_names, cell_covs) {

  # --- 1. DEFINE LOOKUP FIRST ---
  comparison_site_lookup <- clustering_df %>%
    dplyr::select(checklist_id, comparison_site = site)

  # --- 2. FILTER & PREPARE (Updated) ---
  # Use inner_join to implicitly filter train_data to only those checklists 
  # that exist in the clustering_df. 
  train_data_prepped <- train_data %>%
    dplyr::inner_join(comparison_site_lookup, by = "checklist_id")
  
  # --- 3. CLEAN UP COLUMNS ---
  # If 'site' exists (e.g. simulation truth), remove it so we can use comparison_site
  if ("site" %in% names(train_data_prepped)) {
    train_data_prepped <- dplyr::select(train_data_prepped, -site)
  }
  
  train_data_prepped <- train_data_prepped %>%
    dplyr::rename(site = comparison_site) 
  
  # --- 4. CREATE NUMERIC SITE INDICES ---
  M <- nrow(w_matrix)
  
  # Check if w_matrix has row names, otherwise assume 1:M
  if(!is.null(rownames(w_matrix))) {
      site_char_ids <- rownames(w_matrix)
  } else {
      site_char_ids <- as.character(1:M)
  }

  site_id_lookup <- data.frame(
    site_char = site_char_ids,  
    site_numeric = 1:M
  )
  
  # Ensure character matching
  train_data_prepped$site <- as.character(train_data_prepped$site)
  site_id_lookup$site_char <- as.character(site_id_lookup$site_char)
  
  # --- 5. RE-INDEX & CREATE VISIT IDs ---
  train_data_prepped <- train_data_prepped %>%
    dplyr::inner_join(site_id_lookup, by = c("site" = "site_char")) %>%
    dplyr::select(-site) %>%
    dplyr::rename(site = site_numeric) %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(visit_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # --- 6. PIVOT OBSERVATIONS (y) ---
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
  
  # --- 7. PIVOT OBSERVATION COVARIATES ---
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
  
  # --- 8. CREATE UNMARKED FRAME ---
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
                            lower = -Inf, upper = Inf,
                            init_lower = -Inf, init_upper = Inf) {
  
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

    rand_starts <- runif(n_params, min = init_lower, max = init_upper)
    


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