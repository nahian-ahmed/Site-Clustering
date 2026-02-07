###########################
# BayesOptClustGeo helper
# (Updated: GLM Proxy with Validation Set)
###########################

library(rBayesianOptimization) 
library(dplyr)
library(FNN)   # install.packages("FNN") for fast nearest neighbor
library(ClustGeo)

# --- Normalization Helper ---
normalize <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
}

# --- Main Optimization Function ---
bayesianOptimizedClustGeo <- function(
    train_data,          # The "Spatial Input" (2-10 visits)
    validation_data,     # The "Discarded" Data (Singletons + Overflow)
    state_covs, 
    obs_covs,            # Needed for the GLM proxy
    n_init = 15,
    n_iter = 25,
    env_weight = 5       # Weight for environmental distance in ClustGeo
){
    
    # --- 1. PREP DATA ---
    # We only cluster unique locations
    train_locs <- train_data %>%
        dplyr::distinct(locality_id, latitude, longitude, .keep_all = TRUE)
    
    # --- 2. PRE-CALCULATE DISTANCES ---
    # (Calculate once to save time inside the loop)
    df_for_dist <- train_locs
    
    # Normalize features for distance calculation
    for (col in c("latitude", "longitude")) {
        df_for_dist[[col]] <- normalize(df_for_dist[[col]])
    }
    for (col in state_covs) {
        df_for_dist[[col]] <- normalize(df_for_dist[[col]]) * env_weight
    }
    
    env_dist <- dist(df_for_dist[, state_covs])
    geo_dist <- dist(df_for_dist[, c("latitude", "longitude")])
    
    # --- 3. FITNESS FUNCTION ---
    # This runs inside the Bayesian Optimization loop
    clustGeo_fit <- function(rho, kappa) {
        
        # Wrap in tryCatch to handle rare clustering failures
        tryCatch({
            
            # A. Generate Clusters
            # --------------------
            tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = rho)
            
            percent <- kappa / 100.0
            K <- max(2, round(nrow(train_locs) * percent))
            
            # Get Cluster IDs for Training Locations
            train_locs$site_id <- cutree(tree, K)
            
            # B. Assign Validation Data to Clusters
            # -------------------------------------
            # 1. Existing Locations: Join by locality_id
            val_mapped <- validation_data %>%
                left_join(train_locs[, c("locality_id", "site_id")], by = "locality_id")
            
            # 2. New Locations (Singletons): Nearest Neighbor
            val_found <- val_mapped %>% filter(!is.na(site_id))
            val_missing <- val_mapped %>% filter(is.na(site_id))
            
            if (nrow(val_missing) > 0) {
                # Find nearest training location
                nn_idx <- FNN::get.knnx(
                    data = train_locs[, c("latitude", "longitude")],
                    query = val_missing[, c("latitude", "longitude")],
                    k = 1
                )$nn.index
                val_missing$site_id <- train_locs$site_id[nn_idx]
            }
            
            # Combine validation sets
            full_validation <- bind_rows(val_found, val_missing)
            
            # Join site IDs back to full training data
            full_train <- train_data %>%
                left_join(train_locs[, c("locality_id", "site_id")], by = "locality_id")
            
            # C. Calculate SITE-LEVEL Covariates
            # ----------------------------------
            # CRITICAL: We average the environment per cluster.
            # This forces the GLM to evaluate the CLUSTER'S predictive power.
            site_covs <- full_train %>%
                group_by(site_id) %>%
                summarise(across(all_of(state_covs), mean, .names = "site_{.col}"), .groups="drop")
            
            # Attach Site Covs to Train & Validation
            model_train <- full_train %>% left_join(site_covs, by = "site_id")
            model_val <- full_validation %>% left_join(site_covs, by = "site_id")
            
            # D. Fit GLM Proxy
            # ----------------
            # Formula: observed ~ site_environment + observation_covariates
            site_vars <- paste0("site_", state_covs)
            f_str <- paste("species_observed ~", 
                           paste(c(site_vars, obs_covs), collapse = " + "))
            
            # Fit on Training
            m_proxy <- glm(as.formula(f_str), data = model_train, family = binomial())
            
            # Predict on Validation
            preds <- predict(m_proxy, newdata = model_val, type = "response")
            
            # E. Calculate AUC
            # ----------------
            labels <- model_val$species_observed
            
            # Edge case: Validation set has only 0s or only 1s (rare but possible)
            if(length(unique(labels)) < 2) return(list(Score = 0, Pred = 0))
            
            # Simple AUC calculation using PRROC (already loaded in main script) or ROCR
            # Using PRROC since it's in your dependencies
            fg <- preds[labels == 1]
            bg <- preds[labels == 0]
            
            if(length(fg) == 0 || length(bg) == 0) return(list(Score = 0, Pred = 0))
            
            auc_score <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg)$auc
            
            return(list(Score = auc_score, Pred = 0))
            
        }, error = function(e) {
            # Return low score on failure
            return(list(Score = 0, Pred = 0))
        })
    }
    
    # --- 4. RUN OPTIMIZATION ---
    search_bounds <- list(rho = c(0.1, 0.9), kappa = c(10, 90))
    
    search_grid <- data.frame(
        rho = runif(n_init, 0.1, 0.9),
        kappa = runif(n_init, 10, 90)
    )
    
    res <- rBayesianOptimization::BayesianOptimization(
        FUN = clustGeo_fit,
        bounds = search_bounds,
        init_grid_dt = search_grid,
        init_points = 0,
        n_iter = n_iter,
        acq = "ucb",
        verbose = FALSE
    )
    
    return(res)
}