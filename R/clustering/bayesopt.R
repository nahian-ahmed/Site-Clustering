###########################
# BayesOptClustGeo helper
# (Refactored for performance and safety)
#
# November 9, 2025
###########################

library(rBayesianOptimization) # bayesian optimization library
library(cluster)               # silhouette width calculation
library(dplyr)                 # For distinct() and %>%

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
    
    m_dist <- distances::distances(df)
    
    return(m_dist)
}


# --- Main Bayesian Optimization Function ---

bayesianOptimizedClustGeo <- function(
    train_data, 
    state_covs, 
    obs_covs, 
    fit_func,
    n_init = 20,
    n_iter = 10  
){
    
    # --- 1. PRE-CALCULATE ALL STATIC DATA ---
    # This data doesn't change between iterations,
    # so we compute it ONCE.
    
    message("BayesOpt: Pre-processing data...")

    # Create the dataframe for linking full data to unique locations
    train_data_cgul <- train_data
    train_data_cgul$lat_long <- paste0(train_data_cgul$latitude, "_", train_data_cgul$longitude)
    
    # Get the unique locations to be clustered
    uniq_loc_df <- dplyr::distinct(train_data_cgul, lat_long, .keep_all = T)
    
    # Pre-calculate the expensive distance matrix ONCE
    # This matrix is used for the silhouette score
    message("BayesOpt: Pre-calculating pairwise distance matrix...")
    m_dist <- get_pairwise_distances(train_data_cgul, c("latitude","longitude"), state_covs)
    
    # Pre-calculate prevalence if using that fitness function
    prevalence_val <- NA
    if (fit_func == "prevalence") {
        prevalence_val <- mean(train_data_cgul$species_observed)
    }

    # --- 2. DEFINE THE FITNESS FUNCTION (as a Closure) ---
    #
    # By defining `clustGeo_fit` *inside* this function, it gains
    # access to all the pre-calculated data (m_dist, uniq_loc_df, etc.)
    # and function arguments (state_covs, fit_func, etc.)
    # *without* needing any global variables (`<<-`).
    #
    clustGeo_fit <- function(alpha, lambda) {

        percent <- lambda / 100.0
        
        # Run the clustering (this is the part that must be dynamic)
        clustGeo_df_i <- clustGeoSites(
            alpha = alpha, 
            uniq_loc_df, # Use pre-calculated unique locations
            state_covs, 
            obs_covs, 
            num_sites = round(nrow(uniq_loc_df) * percent)
        )
        
        # Link the *full* dataset back to the new cluster IDs
        # We create a simple lookup table from the cluster results
        site_lookup <- clustGeo_df_i[, c("lat_long", "site")]
        
        # This join is much faster than the original for-loop
        joined_data <- dplyr::left_join(train_data_cgul, site_lookup, by = "lat_long")
        current_sites <- joined_data$site
        current_sites[is.na(current_sites)] <- -1 # Handle any NAs
        
        # Calculate Silhouette width using the pre-calculated distance matrix
        m_sil <- silhouette(current_sites, m_dist) 
        silh <- mean(m_sil[, 3], na.rm = TRUE) # Added na.rm for safety
        
        # Return score based on the chosen fitness function
        if (fit_func == "silhouette") {
            result <- list(Score = silh, Pred = 0)
            return(result)
        } 
        else if (fit_func == "prevalence") {
            
            # Combine site IDs and observation status
            temp_df <- data.frame(site = current_sites, 
                                  observed = train_data_cgul$species_observed)
            
            n_sites <- length(unique(temp_df$site))
            
            # Calculate site occupancy rate for this clustering
            site_occ_rate <- temp_df %>%
              group_by(site) %>%
              summarize(has_true = any(observed == TRUE), .groups = 'drop') %>%
              filter(has_true) %>%
              nrow() / n_sites

            sp_objective <- -abs(prevalence_val - site_occ_rate)  
            
            result <- list(Score = sp_objective + silh, Pred = 0)
            return(result)
        }
    } # --- End of nested clustGeo_fit function ---


    # --- 3. DEFINE SEARCH AND RUN OPTIMIZATION ---
    
    # Define the search boundary
    search_bound <- list(alpha = c(0.01, 0.99),
                         lambda = c(10, 90))

    # Define initial search sample
    search_grid <- data.frame(
        alpha = runif(n_init, search_bound$alpha[1], search_bound$alpha[2]),
        lambda = runif(n_init, search_bound$lambda[1], search_bound$lambda[2])
    )
    
    message("BayesOpt: Starting optimization...")
    
    # Bayesian Optimization
    bayesianOptimized <- rBayesianOptimization::BayesianOptimization(
        FUN = clustGeo_fit, # Use the *local* fitness function
        bounds = search_bound, 
        init_grid_dt = search_grid, 
        init_points = 0, # We are using the grid, so no random init points
        n_iter = n_iter, 
        acq = "ucb",
        verbose = TRUE # Good to see progress
    )
    
    message("BayesOpt: Optimization complete.")

    return (list(
        Best_Value = bayesianOptimized$Best_Value,
        Best_Pars = list(
            alpha = bayesianOptimized$Best_Par[1], 
            lambda = bayesianOptimized$Best_Par[2]
        )
    ))
}