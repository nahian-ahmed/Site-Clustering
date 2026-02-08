###########################
# BayesOptClustGeo helper
# (Updated: Returns Final Clusters directly)
###########################

library(rBayesianOptimization) 
library(dplyr)
library(ClustGeo)

bayesianOptimizedClustGeo <- function(
    train_data,          # The "Spatial Input" (2-10 visits)
    validation_data,     # The "Discarded" Data (Singletons + Overflow)
    state_covs, 
    obs_covs,            # Needed for the GLM proxy
    n_init = 20,
    n_iter = 30
){
    
    # --- 1. PREP DATA ---
    train_locs <- train_data %>%
        dplyr::distinct(locality_id, latitude, longitude, .keep_all = TRUE)
    
    # --- 2. PRE-CALCULATE DISTANCES ---
    df_for_dist <- train_locs
    
    # A. Standardize Geography (Z-Score)
    # Matches the scale of environmental covariates
    df_for_dist$latitude  <- as.numeric(scale(df_for_dist$latitude))
    df_for_dist$longitude <- as.numeric(scale(df_for_dist$longitude))
    
    # B. Environment (ALREADY STANDARDIZED)
    # Use Z-scores as-is from input
    
    # Calculate Distances
    env_dist <- dist(df_for_dist[, state_covs])
    geo_dist <- dist(df_for_dist[, c("latitude", "longitude")])
    
    # --- 3. FITNESS FUNCTION ---
    clustGeo_fit <- function(rho, kappa) {
        
        tryCatch({
            # A. Generate Clusters
            tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = rho)
            
            percent <- kappa / 100.0
            K <- max(2, round(nrow(train_locs) * percent))
            
            # Assign Site IDs to Training Locations
            train_locs$site_id <- cutree(tree, K)
            
            # Map Site IDs back to full training data
            full_train <- train_data %>%
                left_join(train_locs[, c("locality_id", "site_id")], by = "locality_id")
            
            # B. Prepare Training Data (Cluster Averages)
            train_site_covs <- full_train %>%
                group_by(site_id) %>%
                summarise(across(all_of(state_covs), mean, .names = "site_{.col}"), .groups="drop")
            
            model_train <- full_train %>% 
                left_join(train_site_covs, by = "site_id")
            
            # C. Prepare Validation Data (Self-Sites)
            model_val <- validation_data %>%
                mutate(across(all_of(state_covs), ~ .x, .names = "site_{.col}"))
            
            # D. Fit GLM Proxy
            site_vars <- paste0("site_", state_covs)
            f_str <- paste("species_observed ~", 
                           paste(c(site_vars, obs_covs), collapse = " + "))
            
            m_proxy <- glm(as.formula(f_str), data = model_train, family = binomial())
            
            # Predict & Score
            preds <- predict(m_proxy, newdata = model_val, type = "response")
            labels <- model_val$species_observed
            
            if(length(unique(labels)) < 2) return(list(Score = 0, Pred = 0))
            
            fg <- preds[labels == 1]
            bg <- preds[labels == 0]
            
            if(length(fg) == 0 || length(bg) == 0) return(list(Score = 0, Pred = 0))
            
            auc_score <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg)$auc
            
            return(list(Score = auc_score, Pred = 0))
            
        }, error = function(e) {
            return(list(Score = 0, Pred = 0))
        })
    }
    
    # --- 4. RUN OPTIMIZATION ---
    search_bounds <- list(rho = c(0.01, 0.99), kappa = c(5, 95))
    search_grid <- data.frame(
        rho = runif(n_init, 0.01, 0.99),
        kappa = runif(n_init, 5, 95)
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
    
    # --- 5. GENERATE FINAL CLUSTERS ---
    # Retrieve best parameters
    best_rho <- res$Best_Par["rho"]
    best_kappa <- res$Best_Par["kappa"]
    
    # Rerun ONE last time to get the object
    final_tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = best_rho)
    final_percent <- best_kappa / 100.0
    final_K <- max(2, round(nrow(train_locs) * final_percent))
    
    train_locs$site_id <- cutree(final_tree, final_K)
    
    # Return everything in a list
    return(list(
        Optimization = res,
        Final_Site_Map = train_locs[, c("locality_id", "site_id")] # Just the mapping
    ))
}