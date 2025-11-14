library(dplyr)
library(rje) # For expit()
library(terra) # Added for raster operations
library(sf) # Added for spatial operations
library(Matrix) # Added for sparse 'w' matrix

source(file.path("R","utils.R"))
source(file.path("R","clustering_helpers.R"))
source(file.path("R","model_helpers.R"))


prepare_train_data <- function (
    state_covs, 
    obs_covs,
    cov_tif,  
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

  # Now, join the original data with the new, numeric covariate data
  train_df <- inner_join(train_df, train_env_df, by = "checklist_id")
  
  # Use standardized var names
  norm_res <- norm_ds(train_df, obs_covs, state_covs) 

  train_df <- norm_res$df
  norm_list <- norm_res$n_l

  train_df$species_observed <- -1
  train_df$occupied_prob <- -1
  train_df$det_prob <- -1
   
  train_df$formatted_date <- train_df$observation_date

  return (list(train_df = train_df, norm_list = norm_list))
}


simulate_train_data <-  function (
    reference_clustering_df,
    site_geoms_sf,
    parameter_set_row, 
    state_cov_names, 
    obs_cov_names,
    cov_tif,
    norm_list
) {
  
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

  test_filename <- paste0(placeholder_spec_name, "_zf_filtered_region_2018.csv")
  test_df_og <- read.delim(
    file.path("checklist_data","species", placeholder_spec_name, test_filename), 
    sep = ",", header = T
  )


  test_df_og <- test_df_og[!is.na(test_df_og$duration_minutes),]
  test_df_og <- test_df_og[
    test_df_og$observation_date >= "2018-05-15" & 
    test_df_og$observation_date <= "2018-07-09",
  ]
  test_df <- test_df_og
  
  # === 2. EXTRACT & NORMALIZE COVARIATES ===
  test_env_df <- extract_state_covs(test_df, cov_tif) # Use new function name
  test_df <- inner_join(test_df, test_env_df, by = "checklist_id")
  
  # Use standardized var names
  norm_res <- norm_ds(test_df, obs_cov_names, state_cov_names, norm_list = norm_list)
  
  test_df <- norm_res$df
  
  # === 3. EXTRACT PARAMETERS ===
  state_par_list <- as.list(parameter_set_row[, c("state_intercept", state_cov_names)])
  obs_par_list <- as.list(parameter_set_row[, c("obs_intercept", obs_cov_names)])
  
  names(state_par_list)[1] <- "intercept"
  names(obs_par_list)[1] <- "intercept"
  

  # === 4. SIMULATE OCCUPANCY & DETECTION (CHECKLIST-LEVEL) ===
  # For test data, we simulate everything at the checklist-level
  
  test_df$occupied_prob <- calculate_weighted_sum(state_par_list, test_df)
  test_df$occupied_prob <- rje::expit(test_df$occupied_prob)
  test_df$occupied <- rbinom(nrow(test_df), 1, test_df$occupied_prob)

  test_df$det_prob <- calculate_weighted_sum(obs_par_list, test_df)
  test_df$det_prob <- rje::expit(test_df$det_prob)
  test_df$detection <- rbinom(nrow(test_df), 1, test_df$det_prob)

  
  # === 5. FINAL OBSERVATION ===
  test_df$species_observed <- test_df$occupied * test_df$detection

  return (test_df)
}