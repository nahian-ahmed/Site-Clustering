library(dplyr)
library(rje) # For expit()

source(file.path("R","utils.R"))
source(file.path("R","clustering_helpers.R"))
source(file.path("R","model_helpers.R"))


prepare_train_data <- function (
    cov_tif, 
    state_covs, 
    obs_covs, 
    placeholder_spec_name = "AMCR"
){

  train_filename <- paste0(placeholder_spec_name, "_zf_filtered_region_2017.csv")
  train_df_og <- read.delim(
    file.path("checklist_data","species", placeholder_spec_name, train_filename), 
    sep = ",", header = T
  )
  
  train_df_og <- train_df_og[!is.na(train_df_og$duration_minutes),]
  train_df_og <- train_df_og[
    train_df_og$observation_date >= "2017-05-15" & 
    train_df_og$observation_date <= "2017-07-09",
  ]
   
  train_df <- train_df_og

  # Use the standardized function name from utils.R
  train_env_df <- extract_state_covs(train_df, cov_tif) 
  train_df <- inner_join(train_df, train_env_df, by = "checklist_id")
  
  # Use standardized var names
  norm_res <- norm_ds(train_df, obs_covs, state_covs) 

  train_df_unnorm <- train_df[,c(obs_covs, state_covs)]

  train_df <- norm_res$df
  norm_list <- norm_res$n_l
  
  for(name in c(obs_covs, state_covs)){
    train_df[,paste("unnorm_",name)] = train_df_unnorm[,name]
  }

  train_df$species_observed <- -1
  train_df$occupied_prob <- -1
  train_df$det_prob <- -1
   
  train_df$formatted_date <- train_df$observation_date

  return (list(train_df = train_df, norm_list = norm_list))
}



simulate_train_data <-  function (
    base_train_df, 
    clustering_method_name, 
    parameter_set_row, 
    state_cov_names, 
    obs_cov_names
) {
  
  # === 1. GET REFERENCE CLUSTERING ===
  # Use the dispatcher to run the *single* reference clustering method
  ref_clustering_list <- get_clusterings(
      method_names = list(clustering_method_name), # Must be a list
      og_data = base_train_df,
      state_covs = state_cov_names,
      obs_covs = obs_cov_names,
      truth_df = NULL # No truth_df when generating the reference
  )
  
  # Get the resulting data frame (it's the first and only item in the list)
  sites_df <- ref_clustering_list[[1]] 
  
  
  # === 2. EXTRACT PARAMETERS ===
  # Directly extract parameters from the provided data frame row
  occ_par_list <- as.list(parameter_set_row[, c("state_intercept", state_cov_names)])
  det_par_list <- as.list(parameter_set_row[, c("obs_intercept", obs_cov_names)])
  
  # Rename intercept keys to match what calculate_weighted_sum expects
  names(occ_par_list)[1] <- "intercept"
  names(det_par_list)[1] <- "intercept"
  
  
  # === 3. ENFORCE CLOSURE ===
  sites_list <- unique(sites_df$site)
  closed_df <- enforceClosure(sites_df, state_cov_names, sites_list)
  
  
  # === 4. SIMULATE OCCUPANCY (SITE-LEVEL) ===
  # Get unique sites *from the closed_df* to ensure mean covariates are used
  sites_df_u <- subset(closed_df, !duplicated(site)) 
  
  sites_df_u$occupied_prob <- calculate_weighted_sum(occ_par_list, sites_df_u)
  sites_df_u$occupied_prob <- rje::expit(sites_df_u$occupied_prob)
  sites_df_u$occupied <- rbinom(nrow(sites_df_u), 1, sites_df_u$occupied_prob)
  
  # Create a simple lookup for site-level values
  site_lookup <- sites_df_u[, c("site", "occupied_prob", "occupied")]
  
  
  # === 5. MAP SITE-LEVEL OCCUPANCY TO CHECKLISTS ===
  # A dplyr join is much faster and cleaner than the R-style for-loop
  res_df <- dplyr::left_join(closed_df, site_lookup, by = "site")
  
  
  # === 6. SIMULATE DETECTION (CHECKLIST-LEVEL) ===
  res_df$det_prob <- calculate_weighted_sum(det_par_list, res_df)
  res_df$det_prob <- rje::expit(res_df$det_prob)
  res_df$detection <- rbinom(nrow(res_df), 1, res_df$det_prob)
  
  
  # === 7. FINAL OBSERVATION ===
  res_df$species_observed <- res_df$occupied * res_df$detection
  
  
  # === 8. UN-NORMALIZE COVARIATES ===
  cov_names <- c(state_cov_names, obs_cov_names)
  res_df <- res_df[,!names(res_df) %in% cov_names] # Drop normalized columns
  
  # Rename 'unnorm_' columns back to original names
  for(name in cov_names){
    res_df[,name] = res_df[,paste("unnorm_",name)]
    # Drop the 'unnorm_' column
    res_df[,paste("unnorm_",name)] <- NULL 
  }
  
  return (res_df)
}



simulate_test_data <- function (
    norm_list, 
    parameter_set_row, 
    state_cov_names, 
    obs_cov_names, 
    cov_tif, # <-- Pass the raster in
    placeholder_spec_name = "AMCR"
){
  
  # === 1. LOAD & FILTER TEST DATA ===
  f.name_test <- paste0(placeholder_spec_name, "/", placeholder_spec_name, "_zf_filtered_region_2018.csv")
  test.df.og <- read.delim(
    file.path("checklist_data","species", placeholder_spec_name, f.name_test), 
    sep = ",", header = T
  )
  
  test.df.og <- test.df.og[!is.na(test.df.og$duration_minutes),]
  test.df.og <- test.df.og[
    test.df.og$observation_date >= "2018-05-15" & 
    test.df.og$observation_date <= "2018-07-09",
  ]
  test.df <- test.df.og
  
  # === 2. EXTRACT & NORMALIZE COVARIATES ===
  test_env.df <- extract_state_covs(test.df, cov_tif) # Use new function name
  test.df <- inner_join(test.df, test_env.df, by = "checklist_id")
  
  # Use standardized var names
  norm_res <- norm_ds(test.df, obs_cov_names, state_cov_names, test=TRUE, norm.list = norm_list)
  
  test.df_unnorm <- test.df[,c(obs_cov_names, state_cov_names)]
  test.df <- norm_res$df

  # Store unnormalized values
  for(name in c(obs_cov_names, state_cov_names)){
    test.df[,paste("unnorm_",name)] = test_df_unnorm[,name]
  }
  
  
  # === 3. EXTRACT PARAMETERS ===
  occ_par_list <- as.list(parameter_set_row[, c("state_intercept", state_cov_names)])
  det_par_list <- as.list(parameter_set_row[, c("obs_intercept", obs_cov_names)])
  
  names(occ_par_list)[1] <- "intercept"
  names(det_par_list)[1] <- "intercept"
  

  # === 4. SIMULATE OCCUPANCY & DETECTION (CHECKLIST-LEVEL) ===
  # For test data, we simulate everything at the checklist-level
  
  test.df$occupied_prob <- calculate_weighted_sum(occ_par_list, test.df)
  test.df$occupied_prob <- rje::expit(test.df$occupied_prob)
  test.df$occupied <- rbinom(nrow(test.df), 1, test.df$occupied_prob)

  test.df$det_prob <- calculate_weighted_sum(det_par_list, test.df)
  test.df$det_prob <- rje::expit(test.df$det_prob)
  test.df$detection <- rbinom(nrow(test.df), 1, test.df$det_prob)

  
  # === 5. FINAL OBSERVATION ===
  test.df$species_observed <- test.df$occupied * test.df$detection

  
  # === 6. UN-NORMALIZE COVARIATES ===
  cov_names <- c(state_cov_names, obs_cov_names)
  test.df <- test.df[,!names(test.df) %in% cov_names]
  
  for(name in cov_names){
    test.df[,name] = test.df[,paste("unnorm_",name)]
    test.df[,paste("unnorm_",name)] <- NULL
  }

  return (test.df)
}