# -----------------------------------------------------------------
# Simulation for occuN model
#
# -----------------------------------------------------------------

############
# 1. Installation
############

install_now = FALSE
if (install_now){
  # Set a CRAN mirror
  options(repos = c(CRAN = "https://cloud.r-project.org/"))

  packages <- c("devtools")
  for(p in packages){
      if (!requireNamespace(p, quietly = FALSE)) install.packages(p)
  }

  # Quietly install forked 'unmarked' package
  suppressMessages(
      devtools::install_github("anonymous97331/unmarked", ref = "main", force = TRUE, quiet = FALSE)
  )

  # --- Load required libraries ---
  library(unmarked)
}

library(dplyr)


source(file.path("R", "utils.R"))
source(file.path("R", "simulation_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))

set.seed(123) # For reproducibility


comparison_method_list <- c(
    "1to10",
    "2to10",
    "2to10-sameObs",
    "1-kmSq",
    "lat-long",
    "rounded-4",
    "SVS",
    "1-per-UL",
    "DBSC",
    "BayesOptClustGeo"
)

comparison_method_list <- c(
    "SVS"
)



# optimizers = "BFGS", "L-BFGS-B", "CG", "Nelder-Mead", "SANN", "nlminb" 
selected_optimizer <- "nlminb"



###########
# 2. LOAD CONFIGURATION
###########
sim_params <- read.delim(file.path("config","simulation_parameters.csv"), sep = ",", header = T)
sim_clusterings <- read.delim(file.path("config","simulation_clusterings.csv"), sep = ",", header = T)


n_simulations <- 25
n_fit_repeats <- 25
n_test_repeats <- 25

n_simulations <- 3
n_fit_repeats <- 3
n_test_repeats <- 3



state_cov_names <- names(sim_params)[2:6]
obs_cov_names <- names(sim_params)[8:12]



# Load the base landscape raster
state_cov_raster <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster) <- state_cov_names

base_train_data <- prepare_train_data(state_cov_names, obs_cov_names, state_cov_raster)
base_train_df <- base_train_data$train_df
norm_list <- base_train_data$norm_list



reference_method_list <- sim_clusterings$method
all_method_names <- unique(c(reference_method_list, comparison_method_list))

print(all_method_names)

cat("--- Pre-computing ALL clusterings... ---\n")
all_clusterings <- get_clusterings(
    method_names = all_method_names,
    og_data = base_train_df,
    state_covs = state_cov_names,
    truth_df = NULL # truth_df is now just another method
)
cat(sprintf("--- Pre-computing complete. Found %d total clusterings. ---\n", length(all_clusterings)))


cat("--- Pre-computing site geometries for ALL clustering methods... ---\n")
all_site_geometries <- list() # Renamed list for clarity
for (method_name in all_method_names) {
    cat(paste("    - Generating geometries for:", method_name, "\n"))
    
    # Get the correct data (handles BayesOpt list structure)
    cluster_data <- all_clusterings[[method_name]]
    if (is.list(cluster_data) && "result_df" %in% names(cluster_data)) {
      cluster_data <- cluster_data$result_df
    }
    
    # Handle cases where a method might have failed (e.g., reference_clustering)
    if (is.null(cluster_data)) {
        cat(paste("    - WARNING: No clustering data found for", method_name, ". Skipping geometry creation.\n"))
        next
    }
    
    all_site_geometries[[method_name]] <- create_site_geometries(
        cluster_data, 
        state_cov_raster
    )
}
cat("--- Geometry pre-computing complete. ---\n")

# Calculate summaries using "m" (meters) since raster/geometries are in Albers
# Change units = "km" if you prefer kilometers/sq km
clustering_summary_df <- summarize_clusterings(
  all_clusterings = all_clusterings,
  all_site_geometries = all_site_geometries,
  units = "m" 
)

# Define output directory (create if it doesn't exist yet)
output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the summary to CSV
write.csv(clustering_summary_df, file.path(output_dir, "clustering_descriptive_stats.csv"), row.names = FALSE)
cat(sprintf("--- Summary metrics saved to %s/clustering_descriptive_stats.csv ---\n", output_dir))


all_dataset_stats <- list()
all_results <- list()


for (cluster_idx in seq_len(nrow(sim_clusterings))) {

  current_clustering_method <- sim_clusterings$method[cluster_idx]
  print(paste("STARTING REFERENCE CLUSTERING:", current_clustering_method))
  
  current_reference_dataframe <- all_clusterings[[current_clustering_method]]

  current_site_geometries <- all_site_geometries[[current_clustering_method]]
  
  cat(paste("    - Pre-calculating Albers projection for:", current_clustering_method, "\n"))
  albers_crs_str <- sf::st_crs(current_site_geometries)$wkt

  # Project the base raster *once*
  cov_tif_albers <- terra::project(state_cov_raster, albers_crs_str, method="bilinear", res = 30)

  # Calculate cell area raster *once*
  area_j_raster <- terra::cellSize(cov_tif_albers, unit="m")
  cat(paste("    - Pre-calculation complete.\n"))

  # Loop over reference parameters (iterating by row index)
  for (param_idx in seq_len(nrow(sim_params))) {

    current_parameter_set <- sim_params[param_idx, ]
    print(paste0("STARTING PARAMETER SET ", param_idx, ":"))
    print(current_parameter_set)
  
    # Loop for each of the 100 simulations (stochastic replicates)
    for (sim_num in 1:n_simulations) {
      
      print(paste("--- Simulation run", sim_num, "of", n_simulations, "---"))

      # === 1. SIMULATE DATA ===
      # NEW CALL 1:
      train_data <- simulate_train_data(
        reference_clustering_df = current_reference_dataframe, 
        site_geoms_sf = current_site_geometries,
        parameter_set_row = current_parameter_set, 
        state_cov_names = state_cov_names, 
        obs_cov_names = obs_cov_names,
        # cov_tif = state_cov_raster, # <-- REMOVED
        norm_list = norm_list,
        cov_tif_albers = cov_tif_albers, # <-- PASS PRE-CALCULATED RASTER
        area_j_raster = area_j_raster    # <-- PASS PRE-CALCULATED RASTER
      )
            
      # NEW CALL 2:
      test_data_full <- simulate_test_data(
          norm_list = norm_list, 
          parameter_set_row = current_parameter_set, 
          state_cov_names = state_cov_names, 
          obs_cov_names = obs_cov_names,
          cov_tif = state_cov_raster
      )

      # === 1.5. CALCULATE DATASET STATS ===
      # Call the new function
      current_dataset_stats <- summarize_datasets(train_data, test_data_full)

      # Add identifiers for this run
      current_dataset_stats$cluster_method <- current_clustering_method
      current_dataset_stats$param_set <- param_idx
      current_dataset_stats$sim_num <- sim_num

      current_dataset_stats <- current_dataset_stats %>%
        dplyr::relocate(cluster_method, param_set, sim_num)


      # Add the row to the storage list
      all_dataset_stats[[length(all_dataset_stats) + 1]] <- current_dataset_stats
      

      # === 2. GET ALL CLUSTERINGS ===
      train_data_for_clustering <- subset(train_data, select = -c(site))
      
      
      # This model_list should be *inside* the sim_num loop
      # to get fresh models for each stochastic replicate.
      model_list <- list()

      # === 3. RUN EXPERIMENTS FOR EACH TEST REPEAT ===
      for (repeat_num in 1:n_test_repeats) {

        # # Spatially subsample the test data for this repeat
        # # (This logic is from your R/utils.R and OLD_code)
        # set.seed(repeat_num) # Ensure this subsample is reproducible
        # hexagons <- dggridR::dgconstruct(spacing = 5, topology = "HEXAGON")
        # test_data_full$cell <- dggridR::dgGEO_to_SEQNUM(hexagons, test_data_full$latitude, test_data_full$longitude)$seqnum
        
        # test.det.df <- test_data_full[test_data_full$species_observed == T,]
        # test.undet.df <- test_data_full[test_data_full$species_observed == F,]
        
        # cell_names <- names(table(test.det.df$cell))
        # det.valid.df <- spatial_subsample(test.det.df, cell_names) # from R/utils.R
        
        # undet.cell_names <- names(table(test.undet.df$cell))
        # undet.valid.df <- spatial_subsample(test.undet.df, undet.cell_names) # from R/utils.R
        
        # if(nrow(undet.valid.df) > nrow(det.valid.df)){
        #     idx <- sample(seq(1:nrow(undet.valid.df)), nrow(det.valid.df))
        #     undet.valid.df <- undet.valid.df[idx,]
        # } else if (nrow(undet.valid.df) < nrow(det.valid.df)){
        #     idx <- sample(seq(1:nrow(det.valid.df)), nrow(undet.valid.df))
        #     det.valid.df <- det.valid.df[idx,]
        # }
        # test.df <- rbind(det.valid.df, undet.valid.df) # This is the final test set for this repeat


        # # === 4. LOOP OVER EACH CLUSTERING METHOD ===
        # # This is the logic adapted from your `OLD_code/simulation_experiments/run_experiments.R`
        # for(method_name in names(all_clusterings)){
        
        #   current_clustering_df <- all_clusterings[[method_name]]
        #   set.seed(1) # Consistent seed for modeling
          
        #   # --- This logic is for the *first* test repeat only ---
        #   if(repeat_num == 1){
            
        #     # Handle BayesOpt special case
        #     if (startsWith(method_name, "BayesOpt")){
        #         best_par_df <- data.frame(parameter = names(groupedSite[[method_name]]$Best_Pars), value = as.numeric(groupedSite[[method_name]]$Best_Pars))
        #         # ... code to write best_par_df to a CSV ...
                
        #         # IMPORTANT: Overwrite the list entry with just the dataframe
        #         current_clustering_df <- all_clusterings[[method_name]]$result_df
        #     }

        #     # Calculate clustering stats (ARI, etc.) against the ground truth
        #     # We use `train_data` as the ground truth here.
        #     if (!(method_name %in% c("1to10", "2to10", "2to10-sameObs", "1-UL"))){
        #         cl_stats <- calcClusteringStats(current_clustering_df, train_data) # from R/analysis_helpers.R
        #     } else {
        #         cl_stats <- list(ari=NA, ami=NA , nid=NA)
        #     }

        #     cl_stats_desc <- calcDescriptiveClusteringStatsWithReference(current_clustering_df, "site", state_cov_names, normalize = FALSE) # from R/analysis_helpers.R
        #     cl_stats <- c(cl_stats, cl_stats_desc)
            
        #     # ... code to write cl_stats to a CSV ...
        #     # ... code to write groupedSite[[method_name]][,c("checklist_id", "site")] to a CSV ...

        #     # Fit the occupancy model
        #     model_list[[method_name]] <- list()
        #     test.formula <- calcOccModel(current_clustering_df, state_cov_names, obs_cov_names) # from R/model_helpers.R

        #     occ_par_list <- test.formula@estimates@estimates$state@estimates 
        #     det_par_list <- test.formula@estimates@estimates$det@estimates

        #     model_list[[method_name]][["occu_parameters"]] <- occ_par_list
        #     model_list[[method_name]][["det_parameters"]] <- det_par_list 
            
        #     # ... code to calculate occ_par_mape, det_par_mape ...
        #     # ... code to create and write model_pars_df to CSV ...
        #   } # --- End of (repeat_num == 1) block ---


        #   # --- This logic runs for *every* test repeat ---
          
        #   # Get the model parameters (fit only on the first repeat)
        #   occ_par_list <- model_list[[method_name]][["occu_parameters"]]
        #   det_par_list <- model_list[[method_name]][["det_parameters"]]

        #   # Calculate predictions on the *current* test.df subsample
        #   test.df$occupied_prob_est <- calculate_weighted_sum(occ_par_list, test.df) # from R/model_helpers.R
        #   test.df$occupied_prob_est <- rje::expit(test.df$occupied_prob_est)
          
        #   test.df$det_prob_est <- calculate_weighted_sum(det_par_list, test.df) # from R/model_helpers.R
        #   test.df$det_prob_est <- rje::expit(test.df$det_prob_est)
          
        #   pred_observ <- unlist(test.df$occupied_prob_est * test.df$det_prob_est)

        #   # ... code to calculate AUC/AUPRC (roc, pr) ...
          
        #   # ... code to calculate occ.mape.i, det.mape.i ...
          
        #   # ... code to create and write predictions_df to CSV ...
          
        #   # ... code to create and write metrics_df to CSV ...

        #   cat(paste("  -", method_name, "AUC:", round(roc$auc, 4), "\n"))

        # } # End loop over comparison methods
      } # End loop over test repeats (repeat_num)


    } # End simulation loop (sim_num)
  } # End parameter loop (param_idx)
} # End clustering loop (cluster_idx)
print("Done")


# After all loops are done (at the end of the script)
final_dataset_stats_df <- dplyr::bind_rows(all_dataset_stats)

# Save the final summary file
write.csv(final_dataset_stats_df, file.path(output_dir, "dataset_descr_stats.csv"), row.names = FALSE)
cat(sprintf("--- Dataset descriptive stats saved to %s/dataset_descr_stats.csv ---\n", output_dir))


# Save the final summary file
write.csv(all_results, file.path(output_dir, "simulation_summary.csv"), row.names = FALSE)
cat(sprintf("--- Simulation summary saved to %s/simulation_summary.csv ---\n", output_dir))

