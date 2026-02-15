###########################
# BayesOptClustGeo helper
# (Final: Spatial Subsampling, No Fallback, Hex Argument)
###########################

library(rBayesianOptimization) 
library(dplyr)
library(ClustGeo)
library(unmarked)
library(terra)
library(sf)

# Ensure dggridR is loaded for spatial subsampling
if(!require(dggridR)) {
  warning("dggridR is required for spatial subsampling but not loaded.")
}

bayesianOptimizedClustGeo <- function(
  train_data,      
  validation_data,   
  state_covs, 
  obs_covs,
  cov_tif_albers,     
  albers_crs,         
  area_raster,        
  buffer_m = 200,     
  hex_m = 100,
  n_iter = 9,
  n_reps = 3,
  stable_reps = 3,
  lower = -Inf,
  upper = Inf,
  init_lower = -Inf,
  init_upper = Inf
){
  
  # --- 1. CONFIGURATION ---
  init_kappa_values <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95)
  search_grid <- data.frame(kappa = init_kappa_values)
  
  # --- 2. PREP DATA (Train) ---
  train_locs <- train_data %>%
    dplyr::distinct(locality_id, latitude, longitude, .keep_all = TRUE)
  
  df_for_dist <- train_locs
  df_for_dist$latitude  <- as.numeric(scale(df_for_dist$latitude))
  df_for_dist$longitude <- as.numeric(scale(df_for_dist$longitude))
  
  env_dist <- dist(df_for_dist[, state_covs])
  geo_dist <- dist(df_for_dist[, c("latitude", "longitude")])
  
  
  # --- 3. PREP DATA (Validation) ---
  cat("  [BayesOpt] Preparing validation data...\n")
  
  # Attach area_j and other spatial structs to validation data ONCE
  val_prep <- prepare_test_spatial_structures(
    test_df = validation_data,
    albers_crs = albers_crs,
    buffer_m = buffer_m,
    cov_raster_albers = cov_tif_albers,
    area_raster = area_raster
  )
  validation_df_ready <- val_prep$test_df
  
  # Diagnostic: Check Validation Truth
  n_obs <- sum(validation_df_ready$species_observed, na.rm=TRUE)
  cat(sprintf("    [Diag] Validation: %d rows, %d detections (%.1f%%)\n", 
              nrow(validation_df_ready), n_obs, (n_obs/nrow(validation_df_ready))*100))
  
  obs_formula <- as.formula(paste("~", paste(obs_covs, collapse = " + ")))
  state_formula <- as.formula(paste("~", paste(state_covs, collapse = " + ")))

  
  # --- 4. FITNESS FUNCTION ---
  clustGeo_fit <- function(kappa) { 
    
    w_matrix <- NULL
    umf <- NULL
    fm <- NULL
    
    tryCatch({
      # --- A. CLUSTERING ---
      tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = 0.5)
      percent <- kappa / 100.0
      K <- max(2, round(nrow(train_locs) * percent))
      train_locs$site <- cutree(tree, K)
      
      clust_df <- train_data %>%
        dplyr::select(-any_of("site")) %>%
        dplyr::inner_join(train_locs[, c("locality_id", "site")], by = "locality_id")
      
      # --- B. GEOMETRIES & W MATRIX ---
      current_geoms <- create_site_geometries(clust_df, cov_tif_albers, buffer_m, "BayesOpt_Temp")
      split_res <- disjoint_site_geometries(current_geoms, clust_df)
      current_geoms <- split_res$geoms
      clust_df <- split_res$data 
      w_matrix <- generate_overlap_matrix(current_geoms, cov_tif_albers)
      
      # --- C. FIT occuN MODEL ---
      full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))[, state_covs, drop = FALSE]
      full_raster_covs[is.na(full_raster_covs)] <- 0
      
      umf <- prepare_occuN_data(train_data, clust_df, w_matrix, obs_covs, full_raster_covs)
      
      # Fit with explicit bounds
      fm <- fit_occuN_model(
        umf, state_formula, obs_formula,
        n_reps = n_reps, stable_reps = stable_reps,
        optimizer = "nlminb",
        lower = lower,
        upper = upper,
        init_lower = init_lower,
        init_upper = init_upper
      )
      
      if (is.null(fm)) {
        return(list(Score = 0, Pred = 0)) 
      }
      
      est_alphas <- coef(fm, 'det')   
      est_betas <- coef(fm, 'state')
      
      # --- D. CLEANUP ---
      rm(w_matrix, umf, fm, current_geoms)
      gc() 
      
      # --- E. SPATIAL SUBSAMPLING VALIDATION ---
      auc_scores <- numeric(25)
      hex_res_km <- hex_m / 1000  # Convert meters to km for dggridR
      
      valid_subsamples <- 0
      
      # Loop 1 to 25 creates 25 reproducible spatial folds
      for (r in 1:25) {
        # 1. Spatial Subsample (using utils.R logic, seeded by r)
        sub_df <- spatial_subsample_dataset(validation_df_ready, hex_res_km, r)
        
        # # 2. Check for Validity (Size & Class Balance)
        # if(nrow(sub_df) < 10 || length(unique(sub_df$species_observed)) < 2) {
        #    auc_scores[r] <- NA 
        #    next
        # }

        # 3. Predict on Subsample
        X_state <- model.matrix(state_formula, data = sub_df)
        X_obs <- model.matrix(obs_formula, data = sub_df)
        
        # Use robust prediction (cap extremes)
        linear_state <- X_state %*% est_betas
        linear_state[linear_state > 20] <- 20 
        linear_state[linear_state < -20] <- -20
        
        pred_psi <- 1 - exp(-(exp(linear_state) * sub_df$area_j))
        pred_det <- plogis(X_obs %*% est_alphas)
        pred_prob <- pred_psi * pred_det
        
        # 4. Calculate AUC
        metrics <- calculate_classification_metrics(pred_prob, sub_df$species_observed)
        
        if (!is.na(metrics$auc)) {
          auc_scores[r] <- metrics$auc
          valid_subsamples <- valid_subsamples + 1
        } else {
          auc_scores[r] <- NA
        }
      }
      
      # Final Score: Mean AUC across valid spatial subsamples
      mean_auc <- mean(auc_scores, na.rm = TRUE)
      
      # If ALL spatial subsamples failed, return 0
      if (is.na(mean_auc) || valid_subsamples == 0) {
         mean_auc <- 0
      }
      
      return(list(Score = mean_auc, Pred = 0))
      
    }, error = function(e) {
      if(exists("w_matrix")) rm(w_matrix)
      if(exists("umf")) rm(umf)
      gc()
      return(list(Score = 0, Pred = 0))
    })
  }
  
  # --- 5. RUN OPTIMIZATION ---
  search_bounds <- list(kappa = c(5, 95))
  
  res <- rBayesianOptimization::BayesianOptimization(
    FUN = clustGeo_fit,
    bounds = search_bounds,
    init_grid_dt = search_grid,
    init_points = 0,
    n_iter = n_iter, 
    acq = "ucb",
    verbose = TRUE
  )
  
  best_kappa <- res$Best_Par["kappa"]
  
  final_tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = 0.5)
  final_percent <- best_kappa / 100.0
  final_K <- max(2, round(nrow(train_locs) * final_percent))
  train_locs$site <- cutree(final_tree, final_K)
  
  return(list(
    Optimization = res,
    Final_Site_Map = train_locs[, c("locality_id", "site")]
  ))
}