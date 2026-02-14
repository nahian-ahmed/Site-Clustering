###########################
# BayesOptClustGeo helper
# (Updated: Direct occuN optimization with Spatial Subsampling)
###########################

library(rBayesianOptimization) 
library(dplyr)
library(ClustGeo)
library(unmarked)
library(terra)
library(sf)

bayesianOptimizedClustGeo <- function(
  train_data,      
  validation_data,   
  state_covs, 
  obs_covs,
  cov_tif_albers,     # NEW: Required for W matrix
  albers_crs,         # NEW: Required for projection
  area_raster,        # NEW: Required for area_j calc
  buffer_m = 200,     # NEW: For geometries
  hex_m = 100,        # NEW: For subsampling
  n_iter = 9,          # NEW: 9 opt rounds + 11 init rounds = 20 total
  n_reps = 3,
  stable_reps = 3
){
  
  # --- 1. CONFIGURATION ---
  # Define the fixed Grid of Kappa values for initialization
  init_kappa_values <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95)
  search_grid <- data.frame(kappa = init_kappa_values)
  
  # --- 2. PREP DATA (Train) ---
  # Unique locations for clustering
  train_locs <- train_data %>%
    dplyr::distinct(locality_id, latitude, longitude, .keep_all = TRUE)
  
  # Pre-calculate distances for ClustGeo (Speedup)
  df_for_dist <- train_locs
  df_for_dist$latitude  <- as.numeric(scale(df_for_dist$latitude))
  df_for_dist$longitude <- as.numeric(scale(df_for_dist$longitude))
  
  env_dist <- dist(df_for_dist[, state_covs])
  geo_dist <- dist(df_for_dist[, c("latitude", "longitude")])
  
  
  # --- 3. PREP DATA (Validation) ---
  # We prepare validation structures ONCE before the loop to save time.
  # This adds 'area_j' and extracts raster covariates for all validation points.
  cat("  [BayesOpt] Preparing validation spatial structures...\n")
  
  val_prep <- prepare_test_spatial_structures(
    test_df = validation_data,
    albers_crs = albers_crs,
    buffer_m = buffer_m,
    cov_raster_albers = cov_tif_albers,
    area_raster = area_raster
  )
  validation_df_ready <- val_prep$test_df
  
  # Define formulas once
  obs_formula <- as.formula(paste("~", paste(obs_covs, collapse = " + ")))
  state_formula <- as.formula(paste("~", paste(state_covs, collapse = " + ")))

  
  # --- 4. FITNESS FUNCTION ---
  clustGeo_fit <- function(kappa) { 
    
    # Initialize cleanup variables
    w_matrix <- NULL
    umf <- NULL
    fm <- NULL
    
    tryCatch({
      # --- A. CLUSTERING ---
      # 1. Generate Clusters (rho fixed at 0.5)
      tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = 0.5)
      
      percent <- kappa / 100.0
      K <- max(2, round(nrow(train_locs) * percent))
      
      # 2. Assign Sites
      train_locs$site <- cutree(tree, K)
      
      # 3. Map back to full training data
      # (Use inner_join to ensure we match rows)
      clust_df <- train_data %>%
        dplyr::select(-any_of("site")) %>% # remove old site col if exists
        dplyr::inner_join(train_locs[, c("locality_id", "site")], by = "locality_id")
      
      
      # --- B. GEOMETRIES & W MATRIX ---
      # 1. Create Polygons
      current_geoms <- create_site_geometries(clust_df, cov_tif_albers, buffer_m, "BayesOpt_Temp")
      
      # 2. Split Disjoint Sites
      split_res <- disjoint_site_geometries(current_geoms, clust_df)
      current_geoms <- split_res$geoms
      clust_df <- split_res$data 
      
      # 3. Generate W Matrix
      w_matrix <- generate_overlap_matrix(current_geoms, cov_tif_albers)
      
      
      # --- C. FIT occuN MODEL ---
      # 1. Prepare Data
      full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))[, state_covs, drop = FALSE]
      full_raster_covs[is.na(full_raster_covs)] <- 0
      
      umf <- prepare_occuN_data(train_data, clust_df, w_matrix, obs_covs, full_raster_covs)
      
      # 2. Fit (Fast settings for optimization: fewer repeats)
      fm <- fit_occuN_model(
        umf, state_formula, obs_formula,
        n_reps = n_reps, stable_reps = stable_reps,
        optimizer = "nlminb"
      )
      
      if (is.null(fm)) return(list(Score = 0, Pred = 0)) # Model failed to converge
      
      # 3. Extract Coefficients
      est_alphas <- coef(fm, 'det')   
      est_betas <- coef(fm, 'state')
      
      
      # --- D. CLEANUP ---
      # CRITICAL: Free memory before validation loop
      rm(w_matrix, umf, fm, current_geoms)
      gc() 
      
      
      # --- E. SPATIAL SUBSAMPLE VALIDATION ---
      auc_scores <- numeric(25)
      
      # Hex resolution in KM for subsampling
      hex_res_km <- hex_m / 1000
      
      for (r in 1:25) {
        # 1. Subsample
        sub_df <- spatial_subsample_dataset(validation_df_ready, hex_res_km, r)
        
        if(nrow(sub_df) < 5) {
           auc_scores[r] <- NA # Skip if too few points
           next
        }

        # 2. Predict PSI (Occupancy)
        # psi = 1 - exp( - exp(Xb) * area )
        X_state <- model.matrix(state_formula, data = sub_df)
        pred_psi <- 1 - exp(-(exp(X_state %*% est_betas) * sub_df$area_j))
        
        # 3. Predict P (Detection)
        X_obs <- model.matrix(obs_formula, data = sub_df)
        pred_det <- plogis(X_obs %*% est_alphas)
        
        # 4. Prob Observed
        pred_prob <- pred_psi * pred_det
        
        # 5. Calculate AUC
        metrics <- calculate_classification_metrics(pred_prob, sub_df$species_observed)
        
        if (!is.na(metrics$auc)) {
          auc_scores[r] <- metrics$auc
        } else {
          auc_scores[r] <- NA
        }
      }
      
      # Final Score: Mean AUC across 25 repeats
      mean_auc <- mean(auc_scores, na.rm = TRUE)
      
      if (is.na(mean_auc)) mean_auc <- 0
      
      return(list(Score = mean_auc, Pred = 0))
      
    }, error = function(e) {
      message(paste("    [Err] Opt failed:", e$message))
      # Ensure cleanup on error
      if(exists("w_matrix")) rm(w_matrix)
      if(exists("umf")) rm(umf)
      gc()
      return(list(Score = 0, Pred = 0))
    })
  }
  
  
  # --- 5. RUN OPTIMIZATION ---
  # Bounds must include the grid range
  search_bounds <- list(kappa = c(5, 95))
  
  res <- rBayesianOptimization::BayesianOptimization(
    FUN = clustGeo_fit,
    bounds = search_bounds,
    init_grid_dt = search_grid,
    init_points = 0,    # 0 because we strictly use init_grid_dt
    n_iter = n_iter,    # 9 additional rounds
    acq = "ucb",
    verbose = TRUE
  )
  
  
  # --- 6. GENERATE FINAL CLUSTERS ---
  # Retrieve best kappa and map back
  best_kappa <- res$Best_Par["kappa"]
  
  # Rerun clustering one last time
  final_tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = 0.5)
  final_percent <- best_kappa / 100.0
  final_K <- max(2, round(nrow(train_locs) * final_percent))
  
  train_locs$site <- cutree(final_tree, final_K)
  
  # Return result structure
  return(list(
    Optimization = res,
    Final_Site_Map = train_locs[, c("locality_id", "site")]
  ))
}