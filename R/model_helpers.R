library(sf)
library(dplyr)
library(terra)


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







create_site_geometries <- function(
  reference_clustering_df, 
  cov_tif
) {
  
  cat("  - (create_site_geometries) Defining CRS...\n")
  albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  wgs84_crs_str <- "+proj=longlat +datum=WGS84"

  cat("  - (create_site_geometries) Calculating buffer distance...\n")
  template_raster_albers <- terra::project(cov_tif, albers_crs_str, res = 30) 
  rast_res_m <- terra::res(template_raster_albers)
  
  print(paste("Cell size (resolution) in meters (x, y):", rast_res_m[1], rast_res_m[2]))
  buffer_dist <- mean(rast_res_m) / 2
  print(paste("Buffer size (half cell) in meters:", buffer_dist))
  cat(sprintf("  - (create_site_geometries) Buffer distance set to: %.2f meters\n", buffer_dist))

  cat("  - (create_site_geometries) Grouping by site to create hulls...\n")
  
  # --- FIX 1: Corrected dplyr/sf logic ---
  site_geoms_albers <- reference_clustering_df %>%
    group_by(site) %>%
    dplyr::summarise(
      .groups = "drop",
      geometry = {
        unique_pts_df <- dplyr::distinct(data.frame(longitude = longitude, latitude = latitude))
        sf_pts <- sf::st_as_sf(
          unique_pts_df, 
          coords = c("longitude", "latitude"), 
          crs = wgs84_crs_str
        )
        hull_wgs84 <- sf::st_convex_hull(sf::st_combine(sf_pts))
        hull_albers <- sf::st_transform(hull_wgs84, crs = albers_crs_str)
        
        # Return the geometry directly, NOT in a list
        sf::st_buffer(hull_albers, dist = buffer_dist)
      }
    )
    # The redundant %>% sf::st_as_sf() pipe is REMOVED
   
  cat(sprintf("  - (create_site_geometries) Generated %d buffered site geometries.\n", nrow(site_geoms_albers)))
  
  # ---
  # --- NEW CODE BLOCK (with FIX 2) ---
  # ---
  
  cat("  - (create_site_geometries) Calculating weight matrix (w) based on cell overlap...\n")

  site_geoms_albers$site_idx <- 1:nrow(site_geoms_albers)
  
  # --- FIX 2: Remove ID = FALSE ---
  w_df <- terra::extract(
    template_raster_albers, 
    site_geoms_albers, 
    weights = TRUE
    # ID = FALSE has been REMOVED. The default (ID = TRUE) is correct.
  )
  
  M <- nrow(site_geoms_albers)
  K <- terra::ncell(template_raster_albers)
  
  # This call will now work because w_df contains the 'site_idx' column
  w_matrix <- Matrix::sparseMatrix(
    i = w_df$site_idx,
    j = w_df$cell,
    x = w_df$weight,
    dims = c(M, K)
  )
  
  attr(site_geoms_albers, "w_matrix") <- w_matrix
  
  cat("  - (create_site_geometries) Weight matrix calculation complete.\n")

  return(site_geoms_albers)
}