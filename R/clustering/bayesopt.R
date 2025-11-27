###########################
# BayesOptClustGeo helper
# (Refactored for performance and safety)
#
# November 9, 2025
###########################

library(rBayesianOptimization) # bayesian optimization library
library(cluster)               # silhouette width calculation
library(dplyr)                 # For distinct() and %>%
library(distances)

# --- Helper Functions ---
# (These are fine as-is)
normalize <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))
}

get_pairwise_distances <- function(df, sp_features, env_features, normalize = TRUE) {
    
    if (normalize){
        for (sp_feature in sp_features){
            df[,sp_feature] <- normalize(df[,sp_feature])
        }
    }
    df <- df[,c(sp_features,env_features)]
    
    m_dist <- distances::distances(as.data.frame(df))
    
    return(m_dist)
}


# --- Main Bayesian Optimization Function ---

# UPDATED SIGNATURE: no longer takes obs_covs
bayesianOptimizedClustGeo <- function(
    train_data, 
    state_covs, 
    fit_func,
    n_init = 30,
    n_iter = 30  
){
    
    # --- 1. PRE-CALCULATE ALL STATIC DATA ---
    message("BayesOpt: Pre-processing data...")

    train_data_cgul <- train_data
    train_data_cgul$lat_long <- paste0(train_data_cgul$latitude, "_", train_data_cgul$longitude)
    
    uniq_loc_df <- dplyr::distinct(train_data_cgul, lat_long, .keep_all = T)
    
    message("BayesOpt: Pre-calculating pairwise distance matrix...")
    m_dist <- get_pairwise_distances(train_data_cgul, c("latitude","longitude"), state_covs)
    
    

    # --- 2. DEFINE THE FITNESS FUNCTION (as a Closure) ---
    clustGeo_fit <- function(rho, kappa) {

        percent <- kappa / 100.0
        

        clustGeo_df_i <- clustGeoSites(
            alpha = rho, 
            uniq_loc_df,
            state_covs = state_covs, 
            num_sites = round(nrow(uniq_loc_df) * percent)
        )
        
        # Link the *full* dataset back to the new cluster IDs
        site_lookup <- clustGeo_df_i[, c("lat_long", "site")]
        
        joined_data <- dplyr::left_join(train_data_cgul, site_lookup, by = "lat_long")
        current_sites <- joined_data$site
        current_sites[is.na(current_sites)] <- -1 
        
        m_sil <- silhouette(current_sites, m_dist) 
        silh <- mean(m_sil[, 3], na.rm = TRUE) 
        
        if (fit_func == "silhouette") {
            result <- list(Score = silh, Pred = 0)
            return(result)
        }
    } # --- End of nested clustGeo_fit function ---


    # --- 3. DEFINE SEARCH AND RUN OPTIMIZATION ---
    search_bound <- list(rho = c(0.01, 0.99),
                         kappa = c(10, 90))

    search_grid <- data.frame(
        rho = runif(n_init, search_bound$rho[1], search_bound$rho[2]),
        kappa = runif(n_init, search_bound$kappa[1], search_bound$kappa[2])
    )
    
    message("BayesOpt: Starting optimization...")
    
    bayesianOptimized <- rBayesianOptimization::BayesianOptimization(
        FUN = clustGeo_fit, 
        bounds = search_bound, 
        init_grid_dt = search_grid, 
        init_points = 0, 
        n_iter = n_iter, 
        acq = "ucb",
        verbose = TRUE 
    )
    
    message("BayesOpt: Optimization complete.")

    return (list(
        Best_Value = bayesianOptimized$Best_Value,
        Best_Pars = list(
            rho = bayesianOptimized$Best_Par[1], 
            kappa = bayesianOptimized$Best_Par[2]
        )
    ))
}