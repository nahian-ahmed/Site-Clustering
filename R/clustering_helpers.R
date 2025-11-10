##############################################
# Clustering Helpers
#
# November 09, 2025
##############################################

# 1. SOURCE MODULAR CLUSTERING LOGIC
source("R/clustering/bayesopt.R")
source("R/clustering/clustgeo.R")
source("R/clustering/dbsc.R")
source("R/clustering/kmsq.R")

# Load required libraries
library(auk)
library(dplyr)
library(plyr) # For round_any

# 2. INTERNAL HELPER FOR CLUSTGEO
# UPDATED CALL: Now correctly passes *only* state_covs
.run_clustgeo_internal <- function(data, alpha, percent, state_covs, obs_covs) {
  
  data_cgul <- data
  data_cgul$lat_long <- paste0(data_cgul$latitude, "_", data_cgul$longitude)
  uniq_loc_df <- dplyr::distinct(data_cgul, lat_long, .keep_all = T)
  
  num_sites <- round(nrow(uniq_loc_df) * percent)
  
  clustGeo_df_i <- clustGeoSites(
    alpha = alpha,
    checklists = uniq_loc_df,
    occ_covs = state_covs, # <-- CORRECTED: Pass state_covs to occ_covs
    num_sites = num_sites
    # det_covs is no longer needed
  )
  
  # Link un-labeled, but lat-long duplicated checklists to the correct site
  data_cgul$site <- -1
  # Use a faster join instead of a for-loop
  site_lookup <- clustGeo_df_i[, c("lat_long", "site")]
  data_cgul <- dplyr::left_join(data_cgul, site_lookup, by = "lat_long", suffix = c("", ".new"))
  
  # Coalesce the new site IDs back into the main 'site' column
  data_cgul$site <- ifelse(!is.na(data_cgul$site.new), data_cgul$site.new, data_cgul$site)
  data_cgul$site.new <- NULL # Clean up
  data_cgul$site[is.na(data_cgul$site)] <- -1 # Ensure no NAs
  
  return(data_cgul)
}


# 3. THE SMART DISPATCHER FUNCTION
run_clustering_method <- function(method_name, og_data, state_covs, obs_covs, truth_df = NULL) {
  
  set.seed(1) # Ensure reproducibility for each method
  
  parts <- strsplit(method_name, "-")[[1]]
  base_method <- parts[1]
  
  # ---
  # A. Parameter-based methods (kmSq, rounded, clustGeo)
  # ---
  
  if (base_method == "kmSq") {
    # Format is "kmSq-[radius_m]", e.g., "kmSq-1000"
    rad_m <- as.double(parts[2])
    result_df <- kmsq.Sites(og_data, rad_m = rad_m)
    
    # Convert config name (kmSq-1000) to canonical name (1-kmSq)
    km_area <- plyr::round_any((rad_m / 1000)^2, 0.125, ceiling)
    canonical_name <- paste0(km_area, "-kmSq")
    
    return(list(name = canonical_name, data = result_df))
    
  } else if (length(parts) == 2 && parts[2] == "kmSq") {
    # *** NEW BLOCK ***
    # Handles "area-style" names like "2-kmSq" from your config
    area_kmSq <- as.numeric(base_method)
    rad_m <- as.integer(sqrt(area_kmSq) * 1000) # Translate area to radius
    
    message(paste("Translating reference method", method_name, "to radius:", rad_m, "m"))
    
    result_df <- kmsq.Sites(og_data, rad_m = rad_m)
    
    # Return with the *original* name as the key
    return(list(name = method_name, data = result_df))
    
  } else if (base_method == "rounded") {
    digits <- as.integer(parts[2])
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
    # DBSC uses *both* state and obs covs in its old implementation, 
    # but the new one seems to only use state_covs inside formatVert.
    # We pass state_covs to its 'occ_covs' param to match dbsc.R
    result_df <- runDBSC(og_data, state_covs, obs_covs)
    return(list(name = "DBSC", data = result_df))
    
  } else if (method_name == "BayesOptClustGeo") {
    bayes_result <- bayesianOptimizedClustGeo(og_data, state_covs, obs_covs, "silhouette")
    alpha <- bayes_result$Best_Pars$alpha
    percent <- bayes_result$Best_Pars$lambda / 100.0
    
    final_df <- .run_clustgeo_internal(og_data, alpha, percent, state_covs, obs_covs)
    
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
        ungroup() 
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
      results[[clustering_result$name]] <- clustering_result$data
      cat(paste("--- Completed. Stored result as:", clustering_result$name, "---\n"))
    }
  }
  
  return(results)
}