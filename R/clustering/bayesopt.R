###########################
# BayesOptClustGeo helper
# (Updated: Bootstrap Validation)
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
  cov_tif_albers,     
  albers_crs,         
  area_raster,        
  buffer_m = 200,     
  n_iter = 9,
  n_reps = 3,
  stable_reps = 3
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
  
  # We still use this helper to attach area_j and raster values efficiently
  val_prep <- prepare_test_spatial_structures(
    test_df = validation_data,
    albers_crs = albers_crs,
    buffer_m = buffer_m,
    cov_raster_albers = cov_tif_albers,
    area_raster = area_raster
  )
  validation_df_ready <- val_prep$test_df
  
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
      
      fm <- fit_occuN_model(
        umf, state_formula, obs_formula,
        n_reps = n_reps, stable_reps = stable_reps,
        optimizer = "nlminb"
      )
      
      if (is.null(fm)) {
        return(list(Score = 0, Pred = 0)) 
      }
      
      est_alphas <- coef(fm, 'det')   
      est_betas <- coef(fm, 'state')
      
      # --- D. CLEANUP ---
      rm(w_matrix, umf, fm, current_geoms)
      gc() 
      
      # --- E. BOOTSTRAP VALIDATION (Replacing Spatial Subsampling) ---
      auc_scores <- numeric(25)
      valid_folds <- 0
      
      # Pre-calculate matrices for the WHOLE validation set once (Speedup)
      X_state_all <- model.matrix(state_formula, data = validation_df_ready)
      X_obs_all   <- model.matrix(obs_formula, data = validation_df_ready)
      
      # Predict for ALL validation points once
      pred_psi_all <- 1 - exp(-(exp(X_state_all %*% est_betas) * validation_df_ready$area_j))
      pred_det_all <- plogis(X_obs_all %*% est_alphas)
      prob_all     <- pred_psi_all * pred_det_all
      truth_all    <- validation_df_ready$species_observed
      
      n_val <- length(prob_all)
      
      for (r in 1:25) {
        # 1. Bootstrap Sample Indices (With Replacement)
        boot_idx <- sample.int(n_val, n_val, replace = TRUE)
        
        # 2. Extract Predictions and Truth for this sample
        prob_boot <- prob_all[boot_idx]
        truth_boot <- truth_all[boot_idx]
        
        # 3. Check for Monotonicity (Must have 0s and 1s)
        if (length(unique(truth_boot)) < 2) {
           auc_scores[r] <- NA
           next
        }
        
        # 4. Calculate AUC
        metrics <- calculate_classification_metrics(prob_boot, truth_boot)
        
        if (!is.na(metrics$auc)) {
          auc_scores[r] <- metrics$auc
          valid_folds <- valid_folds + 1
        } else {
          auc_scores[r] <- NA
        }
      }
      
      # --- F. AGGREGATE ---
      mean_auc <- mean(auc_scores, na.rm = TRUE)
      
      # Final safety check if bootstrapping failed completely (extremely rare)
      if (is.na(mean_auc)) {
         # Fallback to single calculation on full set
         global_metrics <- calculate_classification_metrics(prob_all, truth_all)
         mean_auc <- global_metrics$auc
      }
      
      if (is.na(mean_auc)) mean_auc <- 0
      
      return(list(Score = mean_auc, Pred = 0))
      
    }, error = function(e) {
      message(paste("    [Err] Opt failed:", e$message))
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