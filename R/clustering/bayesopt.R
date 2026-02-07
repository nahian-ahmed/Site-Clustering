###########################
# BayesOptClustGeo helper
# (Updated: GLM Proxy - No FNN, Validation = Self-Sites)
###########################

library(rBayesianOptimization) 
library(dplyr)
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
    # We only cluster unique locations in the training set
    train_locs <- train_data %>%
        dplyr::distinct(locality_id, latitude, longitude, .keep_all = TRUE)
    
    # --- 2. PRE-CALCULATE DISTANCES (Training Only) ---
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
    clustGeo_fit <- function(rho, kappa) {
        
        tryCatch({
            
            # A. Generate Clusters (Training Only)
            # ------------------------------------
            tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = rho)
            
            percent <- kappa / 100.0
            K <- max(2, round(nrow(train_locs) * percent))
            
            # Assign Site IDs to Training Locations
            train_locs$site_id <- cutree(tree, K)
            
            # Map Site IDs back to full training data
            full_train <- train_data %>%
                left_join(train_locs[, c("locality_id", "site_id")], by = "locality_id")
            
            # B. Prepare Training Data for GLM (Cluster Averages)
            # ---------------------------------------------------
            # We average the environment per cluster.
            train_site_covs <- full_train %>%
                group_by(site_id) %>%
                summarise(across(all_of(state_covs), mean, .names = "site_{.col}"), .groups="drop")
            
            model_train <- full_train %>% 
                left_join(train_site_covs, by = "site_id")
            
            # C. Prepare Validation Data for GLM (Self-Sites)
            # -----------------------------------------------
            # Treat every validation point as its own site. 
            # The "Site Environment" is just the point's raw environment.
            # We create columns named 'site_elevation', etc., using the raw values.
            
            model_val <- validation_data %>%
                mutate(across(all_of(state_covs), ~ .x, .names = "site_{.col}"))
            
            # D. Fit GLM Proxy & Predict
            # --------------------------
            # Formula: observed ~ site_environment + observation_covariates
            site_vars <- paste0("site_", state_covs)
            f_str <- paste("species_observed ~", 
                           paste(c(site_vars, obs_covs), collapse = " + "))
            
            # Fit on Training (Clustered View)
            m_proxy <- glm(as.formula(f_str), data = model_train, family = binomial())
            
            # Predict on Validation (Raw View)
            preds <- predict(m_proxy, newdata = model_val, type = "response")
            
            # E. Calculate AUC
            # ----------------
            labels <- model_val$species_observed
            
            # Edge case checks
            if(length(unique(labels)) < 2) return(list(Score = 0, Pred = 0))
            
            fg <- preds[labels == 1]
            bg <- preds[labels == 0]
            
            if(length(fg) == 0 || length(bg) == 0) return(list(Score = 0, Pred = 0))
            
            # Using PRROC (ensure it's loaded in parent script or load here)
            auc_score <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg)$auc
            
            return(list(Score = auc_score, Pred = 0))
            
        }, error = function(e) {
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