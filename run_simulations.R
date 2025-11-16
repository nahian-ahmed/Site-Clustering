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
sim_params <- read.delim(file.path("config", "simulation_parameters.csv"), sep = ",", header = T)
sim_clusterings <- read.delim(file.path("config", "simulation_clusterings.csv"), sep = ",", header = T)


n_simulations <- 25
n_fit_repeats <- 25
n_test_repeats <- 25

n_simulations <- 1
n_fit_repeats <- 1
n_test_repeats <- 1

res_m <- 30

state_cov_names <- names(sim_params)[2:6]
obs_cov_names <- names(sim_params)[8:12]



# Load the base landscape raster
state_cov_raster <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster) <- state_cov_names

base_train_data <- prepare_train_data(state_cov_names, obs_cov_names, state_cov_raster)
base_train_df <- base_train_data$train_df
norm_list <- base_train_data$norm_list

# +++ NEW: PREPARE TEST DATA ONCE +++
base_test_df <- prepare_test_data(state_cov_names, obs_cov_names, state_cov_raster, norm_list)

cat("--- Pre-calculating Albers projection and cell area (ONCE) ---\n")
# Define the Albers CRS string (from R/model_helpers.R)
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# Project ONCE
cov_tif_albers <- terra::project(state_cov_raster, albers_crs_str, method="bilinear", res = res_m)

# Calculate cell size ONCE
area_j_raster <- terra::cellSize(cov_tif_albers, unit="m")
cat("--- Albers projection and cell area complete ---\n")


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

  # Loop over reference parameters (iterating by row index)
  for (param_idx in seq_len(nrow(sim_params))) {

    current_parameter_set <- sim_params[param_idx, ]
    print(paste0("STARTING PARAMETER SET ", param_idx, ":"))
    print(current_parameter_set)

    # +++ ADDED: PRE-CALCULATE N_j_raster and obs_par_list ONCE PER PARAM SET +++
    message("  (main_loop) Pre-calculating N_j_raster...")
    
    # Get state parameters
    state_par_list <- as.list(current_parameter_set[, c("state_intercept", state_cov_names)])
    names(state_par_list)[1] <- "intercept"
    
    # Calculate log_lambda_j_raster
    log_lambda_j_raster <- cov_tif_albers[[1]] * 0 + state_par_list$intercept
    for (cov_name in state_cov_names) {
      if (cov_name %in% names(cov_tif_albers)) {
        log_lambda_j_raster <- log_lambda_j_raster + (cov_tif_albers[[cov_name]] * state_par_list[[cov_name]])
      }
    }
    
    lambda_j_raster <- exp(log_lambda_j_raster)
    
    # Calculate N_j_raster
    N_j_raster <- lambda_j_raster * area_j_raster
    
    # Get obs parameters
    obs_par_list <- as.list(current_parameter_set[, c("obs_intercept", obs_cov_names)])
    names(obs_par_list)[1] <- "intercept"
    
    message("  (main_loop) N_j_raster pre-calculation complete.")
    # +++ END OF ADDED BLOCK +++
  
    # Loop for each of the 100 simulations (stochastic replicates)
    for (sim_num in 1:n_simulations) {
      
      print(paste("--- Simulation run", sim_num, "of", n_simulations, "---"))

      # === 1. SIMULATE DATA ===
      
      # +++ MODIFIED CALL 1 +++
      train_data <- simulate_train_data(
        reference_clustering_df = current_reference_dataframe,
        site_geoms_sf = current_site_geometries,
        # parameter_set_row = current_parameter_set, # Removed
        # state_cov_names = state_cov_names,         # Removed
        obs_cov_names = obs_cov_names,
        obs_par_list = obs_par_list,            # Pass pre-calculated list
        N_j_raster = N_j_raster                # Pass pre-calculated raster
      )
          

      test_data_full <- simulate_test_data(
          base_test_df = base_test_df,
          # parameter_set_row = current_parameter_set, # No longer needed for state
          # state_cov_names = state_cov_names,       # No longer needed for state
          obs_cov_names = obs_cov_names,
          obs_par_list = obs_par_list,            # Pass this instead
          N_j_raster = N_j_raster,                # +++ ADD THIS +++
          albers_crs_str = albers_crs_str
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

        # Spatially subsample the test data for this repeat
        test_df <- spatial_subsample_dataset(test_data_full = test_data_full, spacing = res_m/1000, repeat_num = repeat_num)

        # === 4. LOOP OVER EACH CLUSTERING METHOD ===
        for(method_name in names(all_clusterings)){
        
          current_clustering_df <- all_clusterings[[method_name]]
          print(method_name)

        } # End loop over comparison methods
      } # End loop over test repeats (repeat_num)


    } # End simulation loop (sim_num)
  } # End parameter loop (param_idx)
} # End clustering loop (cluster_idx)
print("Done")


# After all loops are done (at the end of the script)
final_dataset_stats_df <- dplyr::bind_rows(all_dataset_stats)

# Save the final summary file
write.csv(final_dataset_stats_df, file.path(output_dir, "dataset_descriptive_stats.csv"), row.names = FALSE)
cat(sprintf("--- Dataset descriptive stats saved to %s/dataset_descriptive_stats.csv ---\n", output_dir))


# Save the final summary file
write.csv(all_results, file.path(output_dir, "simulation_summary.csv"), row.names = FALSE)
cat(sprintf("--- Simulation summary saved to %s/simulation_summary.csv ---\n", output_dir))