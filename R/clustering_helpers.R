##############################################
# Clustering Helpers
#
# November 09, 2025
##############################################

# 1. SOURCE MODULAR CLUSTERING LOGIC
# (You already do this, which is great!)
source("R/clustering/bayesopt.R")
source("R/clustering/clustgeo.R")
source("R/clustering/dbsc.R")
source("R/clustering/kmsq.R")

# Load required libraries
library(auk)
library(dplyr)

# 2. INTERNAL HELPER FOR CLUSTGEO
# This logic is used by both clustGeo and BayesOptClustGeo,
# so we make it a reusable internal function.
.run_clustgeo_internal <- function(data, alpha, percent, state_covs, obs_covs) {
  
  data_cgul <- data
  data_cgul$lat_long <- paste0(data_cgul$latitude, "_", data_cgul$longitude)
  uniq_loc_df <- dplyr::distinct(data_cgul, lat_long, .keep_all = T)
  
  num_sites <- round(nrow(uniq_loc_df) * percent)
  
  clustGeo_df_i <- clustGeoSites(
    alpha = alpha,
    checklists = uniq_loc_df,
    state_covs = state_covs,
    obs_covs = obs_covs,
    num_sites = num_sites
  )
  
  # Link un-labeled, but lat-long duplicated checklists to the correct site
  data_cgul$site <- -1
  for(j in seq(1:nrow(clustGeo_df_i))){
    data_cgul[data_cgul$lat_long == clustGeo_df_i[j,]$lat_long, ]$site <- clustGeo_df_i[j,]$site
  }
  return(data_cgul)
}


# 3. THE SMART DISPATCHER FUNCTION
# This new function parses the method name and calls the
# correct clustering code. It returns the clustered data
# and the canonical name for the results list.
run_clustering_method <- function(method_name, og_data, state_covs, obs_covs, truth_df = NULL) {
  
  set.seed(1) # Ensure reproducibility for each method
  
  parts <- strsplit(method_name, "-")[[1]]
  base_method <- parts[1]
  
  # ---
  # A. Parameter-based methods (kmSq, rounded, clustGeo)
  # ---
  if (base_method == "kmSq") {
    rad_m <- as.double(parts[2])
    result_df <- kmsq.Sites(og_data, rad_m = rad_m)
    
    # Convert config name (kmSq-1000) to canonical name (1-kmSq)
    # to match your old analysis scripts.
    km_area <- round_any((rad_m / 1000)^2, 0.125, ceiling)
    canonical_name <- paste0(km_area, "-kmSq")
    
    return(list(name = canonical_name, data = result_df))
    
  } else if (base_method == "rounded") {
    digits <- as.integer(parts[2])
    # Assuming you have round_lat_long in R/utils.R
    rounded_data <- round_lat_long(og_data, rounding_degree = digits) 
    result_df <- auk::filter_repeat_visits(
      rounded_data,
      min_obs = 1, max_obs = 1000000, annual_closure = TRUE,
      date_var = "formatted_date", site_vars = c("rounded_locality_id")
    )
    return(list(name = method_name, data = result_df))
    
  } else if (base_method == "clustGeo") {
    alpha <- as.double(parts[2]) / 100.0
    percent <- as.double(parts[3]) / 100.0
    
    result_df <- .run_clustgeo_internal(og_data, alpha, percent, state_covs, obs_covs)
    return(list(name = method_name, data = result_df))
  
  # ---
  # B. Complex, self-contained methods (DBSC, BayesOpt)
  # ---
  } else if (method_name == "DBSC") {
    result_df <- runDBSC(og_data, state_covs, obs_covs)
    return(list(name = "DBSC", data = result_df))
    
  } else if (method_name == "BayesOptClustGeo") {
    # Run Bayesian Optimization to find parameters
    bayes_result <- bayesianOptimizedClustGeo(og_data, state_covs, obs_covs, "silhouette")
    alpha <- bayes_result$Best_Pars$alpha
    percent <- bayes_result$Best_Pars$lambda / 100.0
    
    # Run clustering with the optimized parameters
    final_df <- .run_clustgeo_internal(og_data, alpha, percent, state_covs, obs_covs)
    
    # Special return format for this method, as it also returns the parameters
    return(list(
      name = "BayesOptClustGeo",
      data = list(result_df = final_df, Best_Pars = bayes_result$Best_Pars)
    ))
  
  # ---
  # C. Simple, flag-based methods (auk filters, svs, etc.)
  # ---
  } else {
    auk_params <- switch(method_name,
      "one_to_10" = list(min_obs = 1, max_obs = 10, site_vars = c("locality_id")),
      "two_to_10" = list(min_obs = 2, max_obs = 10, site_vars = c("locality_id")),
      "two_to_10_sameObs" = list(min_obs = 2, max_obs = 10, site_vars = c("locality_id", "observer_id")),
      "lat_long" = list(min_obs = 1, max_obs = 1000000, site_vars = c("locality_id")),
      NULL # Default
    )
    
    if (!is.null(auk_params)) {
      result_df <- auk::filter_repeat_visits(
        og_data,
        min_obs = auk_params$min_obs,
        max_obs = auk_params$max_obs,
        annual_closure = TRUE,
        date_var = "formatted_date",
        site_vars = auk_params$site_vars
      )
      return(list(name = method_name, data = result_df))
    }
    
    # Handle non-auk simple methods
    if (method_name == "svs") {
      result_df <- og_data
      result_df$site <- result_df$checklist_id
      return(list(name = "SVS", data = result_df))
    }
    
    if (method_name == "one_UL") {
      df_1_UL_t <- auk::filter_repeat_visits(
        og_data,
        min_obs = 1, max_obs = 1000000, annual_closure = TRUE,
        date_var = "formatted_date", site_vars = c("locality_id")
      )
      result_df <- df_1_UL_t %>% 
        group_by(site) %>% 
        filter(row_number() == 1) %>%
        ungroup() # Always good to ungroup after
      return(list(name = "1-UL", data = result_df))
    }
    
    if (method_name == "reference_clustering" && !is.null(truth_df)) {
      return(list(name = "reference-clustering", data = truth_df))
    }
  }
  
  # If we got here, the method wasn't found
  warning(paste("No implementation found for method:", method_name))
  return(NULL)
}


# 4. THE NEW, CLEAN getClusterings FUNCTION
# This replaces the old genExp + getClusterings
# It takes the simple list of names from your run_*.R file
# and uses the dispatcher to get the results.

get_clusterings <- function(method_names, og_data, state_covs, obs_covs, truth_df = data.frame()) {
  
  results <- list()
  
  for (method_config_name in method_names) {
    
    cat(paste("--- Running clustering for:", method_config_name, "---\n"))
    
    clustering_result <- run_clustering_method(
      method_config_name,
      og_data,
      state_covs,
      obs_covs,
      truth_df
    )
    
    if (!is.null(clustering_result)) {
      # Use the canonical name returned by the dispatcher as the key
      results[[clustering_result$name]] <- clustering_result$data
      cat(paste("--- Completed. Stored result as:", clustering_result$name, "---\n"))
    }
  }
  
  return(results)
}