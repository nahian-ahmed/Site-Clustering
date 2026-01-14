####################
# Helper functions

# November 9, 2025
####################

library(terra) # geospatial operations
library(dggridR)



# Normalize state covs
norm_state_covs <- function (state_cov_raster_raw){

  # 1. Calculate Min/Max for each layer, explicitly ignoring NAs
  #    'na.rm = TRUE' is critical here.
  r_min <- terra::minmax(state_cov_raster_raw, compute = TRUE)[1,]
  r_max <- terra::minmax(state_cov_raster_raw, compute = TRUE)[2,]

  # 2. Initialize the output raster
  state_cov_raster <- state_cov_raster_raw

  # 3. Loop through bands to normalize safely
  #    This ensures we handle each layer's specific range
  for (i in 1:terra::nlyr(state_cov_raster)) {
    
    # Get range for this specific band
    band_min <- r_min[i]
    band_max <- r_max[i]
    band_range <- band_max - band_min
    
    # Safety check: If band is constant, avoid division by zero
    if (band_range == 0) {
      # If constant, set to 0 (or 0.5, or keep raw if range is 0-1)
      state_cov_raster[[i]] <- 0 
    } else {
      # Apply (x - min) / (max - min)
      # terra automatically preserves NAs during this arithmetic
      state_cov_raster[[i]] <- (state_cov_raster_raw[[i]] - band_min) / band_range
    }
  }
  return(state_cov_raster)
}

#######
# normalize dataset
#######
norm_ds <- function(df, obs_covs, state_covs, norm_list = list()){

  if(length(norm_list) == 0){
    for(name in c(obs_covs, state_covs)){
      # calc mean/var for each cov, if training
      ma <- max(df[[name]])
      mi <- min(df[[name]])
      norm_list[[name]] <- c(ma, mi)
    }
  }
  
  
  # Get column names that are BOTH in the norm_list AND in the dataframe
  cols_to_normalize <- intersect(names(norm_list), names(df))

  # xi - min(x)/(max(x) - min(x))
  # Only loop over the columns that actually exist in df
  for(cov in cols_to_normalize){
    df[[cov]] <- (df[[cov]] - norm_list[[cov]][[2]])/(norm_list[[cov]][[1]] - norm_list[[cov]][[2]])
  }
    
  return(list(df=df, n_l=norm_list))
  
}


standardize_state_covs <- function(state_cov_raster){
  
  # 1. Calculate Global Mean and SD for all layers
  #    returns a dataframe with rows=layers, cols=mean,sd
  r_stats <- terra::global(state_cov_raster, c("mean", "sd"), na.rm = TRUE)
  
  # 2. Initialize params list
  params <- list()
  
  # 3. Store parameters in the format expected by standardize_ds
  #    Structure: list(cov_name = c(mean = x, sd = y))
  for(i in 1:nrow(r_stats)){
    layer_name <- rownames(r_stats)[i]
    params[[layer_name]] <- c(mean = r_stats[i, "mean"], sd = r_stats[i, "sd"])
  }
  
  # 4. Standardize the Raster
  #    (x - mean) / sd
  #    terra handles recycling of the stats vector automatically if dimensions match
  std_raster <- (state_cov_raster - t(r_stats$mean)) / t(r_stats$sd)
  
  return(list(raster = std_raster, params = params))
}


#######
# Standardize dataset (Z-score standardization)
#######
standardize_ds <- function(df, obs_covs, state_covs, standardization_params = list()){
  
  # Combine all columns we intend to scale
  cols_to_scale <- intersect(c(obs_covs, state_covs), names(df))
  
  for(cov in cols_to_scale){
    
    # CASE A: Parameter exists (e.g., from Raster or Training set) -> Apply it
    if(cov %in% names(standardization_params)){
      stats <- standardization_params[[cov]]
      mu <- stats['mean']
      sigma <- stats['sd']
      
      df[[cov]] <- (df[[cov]] - mu) / sigma
      
    } 
    # CASE B: Parameter missing -> Calculate from Data & Store
    else {
      mu <- mean(df[[cov]], na.rm = TRUE)
      sigma <- sd(df[[cov]], na.rm = TRUE)
      
      # Store for return
      standardization_params[[cov]] <- c(mean = mu, sd = sigma)
      
      df[[cov]] <- (df[[cov]] - mu) / sigma
    }
  }
  
  return(list(df = df, standardization_params = standardization_params))
}

#######
# extract environmental features
# at checklist locations
#######
extract_state_covs <- function(df, cov_tif, x = "longitude", y = "latitude", crs_in = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") {
  
  # Create Vector from DF
  df_pts <- vect(df, geom = c(x, y), crs = crs_in)
  
  # Explicitly project to raster CRS if they don't match
  # This suppresses the warning and ensures accuracy
  if (!same.crs(df_pts, cov_tif)) {
    df_pts <- project(df_pts, crs(cov_tif))
  }
  
  # Extract environmental features
  env_vars_df <- data.frame(
    checklist_id = df$checklist_id,
    terra::extract(cov_tif, df_pts)
  )
  return(env_vars_df)
}



#######
# spatial subsampling as defined by: 
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/ddi.13271
#######
spatial_subsample <- function(df, cell_names){
  valid_df <- data.frame()
  i <- 0
  for(freq in table(df$cell)){ 
    i <- i + 1
    if(freq > 1){
      checklists <- df[df$cell == cell_names[i],]
      sample <- checklists[sample(nrow(checklists), 1), ]
    } else {
      sample <- df[df$cell == cell_names[i],]
    }
    valid_df <- rbind(valid_df, sample)
  }  
  return(valid_df)
}


spatial_subsample_dataset <- function(test_data_full, spacing, repeat_num){
  set.seed(repeat_num) # Ensure this subsample is reproducible
  
  hexagons <- dggridR::dgconstruct(spacing = spacing, topology = "HEXAGON")
  test_data_full$cell <- dggridR::dgGEO_to_SEQNUM(hexagons, test_data_full$latitude, test_data_full$longitude)$seqnum
  
  test_det_df <- test_data_full[test_data_full$species_observed == T,]
  test_nondet_df <- test_data_full[test_data_full$species_observed == F,]
  
  cell_names <- names(table(test_det_df$cell))
  det_valid_df <- spatial_subsample(test_det_df, cell_names)
  
  nondet_cell_names <- names(table(test_nondet_df$cell))
  nondet_valid_df <- spatial_subsample(test_nondet_df, nondet_cell_names)
  
  if(nrow(nondet_valid_df) > nrow(det_valid_df)){
      idx <- sample(seq_len(nrow(nondet_valid_df)), nrow(det_valid_df))
      nondet_valid_df <- nondet_valid_df[idx,]
  } else if (nrow(nondet_valid_df) < nrow(det_valid_df)){
      idx <- sample(seq_len(nrow(det_valid_df)), nrow(nondet_valid_df))
      det_valid_df <- det_valid_df[idx,]
  }
  test_df <- rbind(det_valid_df, nondet_valid_df)

  return(test_df)

}


#########
# Rounding Lat/Long
#########
round_lat_long <- function(df, rounding_degree){
  df$rounded_lat <- round(df$latitude, digits = rounding_degree)
  df$rounded_long <- round(df$longitude, digits = rounding_degree)
  df$rounded_locality_id <- paste(as.character(df$rounded_long), as.character(df$rounded_lat), sep = "_")
  
  # Remove the temporary columns by setting them to NULL
  df$rounded_lat <- NULL
  df$rounded_long <- NULL
  
  return(df)
}



# Helper to reconstruct SLIC polygons from the raster
reconstruct_slic_polygons <- function(method_name, train_df, cov_raster) {
  # Parse params from string "SLIC-eta-zeta"
  parts <- strsplit(method_name, "-")[[1]]
  eta <- as.numeric(parts[2])
  zeta <- as.numeric(parts[3])
  
  # Prepare checklist data
  checklists_sf <- sf::st_as_sf(train_df, coords = c("longitude", "latitude"), crs = 4326)
  checklists_sf <- sf::st_transform(checklists_sf, terra::crs(cov_raster))
  
  # 1. Generate Seeds
  seeds <- get_slic_seeds(cov_raster, zeta, checklists_sf, buffer_dist_m = 5000)
  
  if (nrow(seeds) == 0) return(NULL)
  
  # 2. Re-run Clustering to get pixel map
  pixel_df <- perform_slic_clustering(cov_raster, seeds, eta, zeta)
  
  # 3. Convert to Polygons (OPTIMIZED)
  r_temp <- terra::rast(cov_raster, nlyrs=1)
  cell_ids <- terra::cellFromXY(r_temp, as.matrix(pixel_df[, c("x", "y")]))
  r_temp[cell_ids] <- pixel_df$site
  
  # --- FILTERING STEP ---
  # Only keep sites that actually contain a checklist
  occupied_sites <- terra::extract(r_temp, checklists_sf, ID = FALSE)
  valid_site_ids <- unique(na.omit(occupied_sites[,1]))
  
  # Set all other cells to NA so they don't become polygons
  r_temp[!r_temp %in% valid_site_ids] <- NA
  # ----------------------
  
  # Polygonize (Only processes occupied regions now)
  polys <- terra::as.polygons(r_temp, dissolve=TRUE)
  
  if (length(polys) == 0) return(NULL)
  
  sf_polys <- sf::st_as_sf(polys)
  colnames(sf_polys)[1] <- "site"
  sf::st_crs(sf_polys) <- terra::crs(cov_raster)
  
  return(sf_polys)
}