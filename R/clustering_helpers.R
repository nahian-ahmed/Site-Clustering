##############################################
# Clustering Helpers
#
# Updated: January 08, 2026
##############################################

# 1. SOURCE MODULAR CLUSTERING LOGIC
source("R/clustering/kmsq.R")
source("R/clustering/dbsc.R")
source("R/clustering/clustgeo.R")
source("R/clustering/bayesopt.R")
source("R/clustering/slic.R") # Added SLIC


# Load required libraries
library(auk)
library(dplyr)


# 2. INTERNAL HELPER FOR CLUSTGEO
# UPDATED CALL: Now correctly passes *only* state_covs
.run_clustgeo_internal <- function(data, rho, percent, state_covs) {
  
  data_cgul <- data
  data_cgul$lat_long <- paste0(data_cgul$latitude, "_", data_cgul$longitude)
  uniq_loc_df <- dplyr::distinct(data_cgul, lat_long, .keep_all = T)
  
  num_sites <- round(nrow(uniq_loc_df) * percent)
  
  clustGeo_df_i <- clustGeoSites(
    alpha = rho,
    checklists = uniq_loc_df,
    state_covs = state_covs,
    num_sites = num_sites
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


# 3. DISPATCHER FUNCTION
# Updated to accept cov_raster for SLIC methods
run_clustering_method <- function(method_name, og_data, state_covs, truth_df = NULL, cov_raster = NULL) {
  
  set.seed(1) # Ensure reproducibility for each method
  
  parts <- strsplit(method_name, "-")[[1]]
  base_method <- parts[1]
  
  # ---
  # A. Parameter-based methods (kmSq, rounded, clustGeo, SLIC)
  # ---
  
  if (base_method == "kmSq") {
    # Format is "kmSq-[radius_m]", e.g., "kmSq-1000"
    rad_m <- as.double(parts[2])
    result_df <- kmsq_sites(og_data, rad_m = rad_m)
    
    # Convert config name (kmSq-1000) to canonical name (1-kmSq)
    km_area <- plyr::round_any((rad_m / 1000)^2, 0.125, ceiling)
    canonical_name <- paste0(km_area, "-kmSq")
    
    return(list(name = canonical_name, data = result_df))
    
  } else if (length(parts) == 2 && parts[2] == "kmSq") {
    # Handles "area-style" names like "2-kmSq"
    area_kmSq <- as.numeric(base_method)
    rad_m <- as.integer(sqrt(area_kmSq) * 1000) # Translate area to radius
    
    message(paste("Translating reference method", method_name, "to radius:", rad_m, "m"))
    
    result_df <- kmsq_sites(og_data, rad_m = rad_m)
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
    rho <- as.double(parts[2]) / 100.0
    percent <- as.double(parts[3]) / 100.0
    
    result_df <- .run_clustgeo_internal(og_data, rho, percent, state_covs)
    return(list(name = method_name, data = result_df))
  
  } else if (base_method == "SLIC") {
    # Format: SLIC-[eta]-[zeta]
    if (is.null(cov_raster)) stop("SLIC method requires 'cov_raster' argument.")
    
    eta <- as.numeric(parts[2])
    zeta <- as.numeric(parts[3])
    
    result_df <- slicSites(og_data, state_covs, cov_raster, eta, zeta)
    return(list(name = method_name, data = result_df))

  # ---
  # B. Complex, self-contained methods (DBSC, BayesOpt)
  # ---
  } else if (method_name == "DBSC") {
    result_df <- runDBSC(og_data, state_covs)
    return(list(name = "DBSC", data = result_df))
    
  } else if (method_name == "BayesOptClustGeo") {
    bayes_result <- bayesianOptimizedClustGeo(og_data, state_covs, "silhouette")
    rho <- bayes_result$Best_Pars$rho
    percent <- bayes_result$Best_Pars$kappa / 100.0
    
    final_df <- .run_clustgeo_internal(og_data, rho, percent, state_covs)
    
    return(list(
      name = "BayesOptClustGeo",
      data = list(result_df = final_df, Best_Pars = bayes_result$Best_Pars)
    ))
  
  # ---
  # C. Simple, flag-based methods (auk filters, svs, etc.)
  # ---
  } else {
    auk_params <- switch(method_name,
      "one_to_10" = list(min_obs = 1, max_obs = 10, site_vars = c("locality_id"), canonical = "1to10"),
      "1to10" = list(min_obs = 1, max_obs = 10, site_vars = c("locality_id"), canonical = "1to10"),
      
      "two_to_10" = list(min_obs = 2, max_obs = 10, site_vars = c("locality_id"), canonical = "2to10"),
      "2to10" = list(min_obs = 2, max_obs = 10, site_vars = c("locality_id"), canonical = "2to10"),
      
      "two_to_10_sameObs" = list(min_obs = 2, max_obs = 10, site_vars = c("locality_id", "observer_id"), canonical = "2to10-sameObs"),
      "2to10-sameObs" = list(min_obs = 2, max_obs = 10, site_vars = c("locality_id", "observer_id"), canonical = "2to10-sameObs"),
      
      "lat_long" = list(min_obs = 1, max_obs = 1000000, site_vars = c("locality_id"), canonical = "lat-long"),
      "lat-long" = list(min_obs = 1, max_obs = 1000000, site_vars = c("locality_id"), canonical = "lat-long"),
      
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
      return(list(name = auk_params$canonical, data = result_df))
    }
    
    if (method_name == "svs" || method_name == "SVS") {
      result_df <- og_data
      result_df$site <- result_df$checklist_id
      return(list(name = "SVS", data = result_df)) 
    }
    
    if (method_name == "one_UL" || method_name == "1-per-UL") {
      df_1_UL_t <- auk::filter_repeat_visits(
        og_data,
        min_obs = 1, max_obs = 1000000, annual_closure = TRUE,
        date_var = "formatted_date", site_vars = c("locality_id")
      )
      result_df <- df_1_UL_t %>% 
        group_by(site) %>% 
        filter(row_number() == 1) %>%
        ungroup() 
      return(list(name = "1-per-UL", data = result_df)) 
    }
    
    if (method_name == "reference_clustering" && !is.null(truth_df)) {
      return(list(name = "reference-clustering", data = truth_df))
    }
  }
  
  warning(paste("No implementation found for method:", method_name))
  return(NULL)
}

# 4. THE NEW, CLEAN getClusterings FUNCTION
# Updated to accept cov_raster
get_clusterings <- function(method_names, og_data, state_covs, truth_df = data.frame(), cov_raster = NULL) {
  
  results <- list()
  
  for (method_config_name in method_names) {
    
    cat(paste("--- Running clustering for:", method_config_name, "---\n"))
    
    clustering_result <- run_clustering_method(
      method_config_name,
      og_data,
      state_covs,
      truth_df,
      cov_raster
    )
    
    if (!is.null(clustering_result)) {
      results[[clustering_result$name]] <- clustering_result$data
      cat(paste("--- Completed. Stored result as:", clustering_result$name, "---\n"))
    }
  }
  
  return(results)
}