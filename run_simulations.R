# -----------------------------------------------
# 1. SETUP
# -----------------------------------------------
print("--- 1. Setting up ---")
# Load libraries
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# # --- Install/Load devtools ---
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }

# # --- Install forked 'unmarked' package ---
# # This ensures your occuN version is available
# print("Installing/updating unmarked package with occuN...")
# suppressMessages(
#   devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE, quiet = FALSE)
# )

# --- Load all required libraries ---
library(unmarked) # Your custom occuN version
# library(tidyverse)
library(terra)
library(stringr)
# Add any other libraries your helpers need (e.g., ClustGeo, sf, igraph)
# library(sf)
# library(igraph)
# library(ClustGeo)

# --- Source all your new helper functions ---
# We assume these files now contain the refactored functions
# from your OLD_code/ and occuN_demo.R
print("Loading helper functions...")
source("R/utils.R")             # For: pivot_checklists, norm_ds, etc.
source("R/simulation_helpers.R") # For: prep_train_data, create_w_matrix, generate_simulated_observations
source("R/clustering_helpers.R") # For: apply_clustering (which wraps kmsq, DBSC, etc.)
source("R/model_helpers.R")      # For: create_unmarked_frame, fit_occuN, calculate_metrics

# -----------------------------------------------
# 2. LOAD CONFIGURATION
# -----------------------------------------------
print("--- 2. Loading Configuration ---")
# --- Load the 9 parameter sets ---
# Note: Your CSV seems to be tab-separated, so we use read_tsv
sim_params <- read_tsv("config/simulation_params.csv")

# --- Load the 9 clustering methods ---
sim_clusterings <- read_csv("config/simulation_clusterings.csv")

# --- Create the 81 "species" scenarios ---
scenarios <- cross_join(sim_params, sim_clusterings)

# --- Define simulation loops ---
N_SIMULATIONS <- 100 # "multiple simulations (100)"
N_REPEATS <- 30     # "repeats (30)" -> Used as random starts for occuN
N_TEST_REPEATS <- 25  # "test set repeats (25)"

# --- Get covariate names from the config file ---
# Occupancy (state) covariates (e.g., elevation, TCB, ...)
occ_cov_names <- names(sim_params)[2:6]
# Detection covariates (e.g., day_of_year, duration_minutes, ...)
det_cov_names <- names(sim_params)[8:12]

# -----------------------------------------------
# 3. LOAD BASE DATA (Run Once)
# -----------------------------------------------
print("--- 3. Loading Base Data (Raster & Checklists) ---")

# --- Load the base landscape raster ---
# This raster provides the cell-level covariates for the *entire* study area
# We assume this path is correct relative to the project root.
OR_normalized <- terra::rast("occupancy_feature_raster/occupancy_features_normalized.tif")
names(OR_normalized) <- c("elevation", "TCB", "TCG", "TCW", "TCA") # Ensure names match config
cellCovs_df <- as.data.frame(OR_normalized, cells = TRUE, na.rm = FALSE)
# Rename the 'cell' column to 'cell_id' for clarity
names(cellCovs_df)[1] <- "cell_id"

# --- Load the base checklist template ---
# This uses the `prep_train_data` function (which you'll move from
# OLD_code/simulation_experiments/simulate_data.R to R/simulation_helpers.R).
# This provides the set of *all possible observation locations* (checklists).
checklist_template_data <- prep_train_data(
  occ_covs = occ_cov_names, 
  det_covs = det_cov_names
)
checklist_template <- checklist_template_data$train.df

print(paste("Base data loaded:", nrow(cellCovs_df), "cells and", nrow(checklist_template), "template checklists."))

# -----------------------------------------------
# 4. INITIALIZE RESULTS
# -----------------------------------------------
all_results <- list()

# -----------------------------------------------
# 5. RUN MASTER LOOP
# -----------------------------------------------
print(paste("--- 5. Starting Master Loop for", nrow(scenarios), "Scenarios ---"))




# # Loop over each of the 81 scenarios (species)
# for (i in 1:nrow(scenarios)) {
  
#   scenario_row <- scenarios[i, ]
  
#   # --- 5.1. Get Scenario-Specific Info ---
#   clustering_method <- scenario_row$method
  
#   # Get the 1xN vector of occupancy parameters
#   occ_params <- unlist(scenario_row[, c("occ_intercept", occ_cov_names)])
  
#   # Get the 1xN vector of detection parameters
#   det_params <- unlist(scenario_row[, c("det_intercept", det_cov_names)])
  
#   # Create the model formula for occuN
#   state_formula_chr <- paste(occ_cov_names, collapse = " + ")
#   det_formula_chr <- paste(det_cov_names, collapse = " + ")
#   full_formula <- as.formula(paste("~", det_formula_chr, "~", state_formula_chr))
  
#   print(paste0("--- Starting Scenario ", i, "/", nrow(scenarios), ": Method='", clustering_method, "' ---"))
  
  
#   # --- 5.2. PREPROCESS / CLUSTER (Done ONCE per scenario) ---
#   # This function applies the correct clustering (kmsq, DBSC, etc.)
#   # to the checklist template to assign a 'site' ID to each checklist.
#   # Logic comes from OLD_code/helper/kmsq.R, DBSCHelper.R, clustGeoHelper.R
#   print(paste("  [Scenario", i, "] 1. Applying clustering..."))
#   clustered_checklists <- apply_clustering(checklist_template, clustering_method)
  
  
#   # --- 5.3. CREATE W MATRIX (Done ONCE per scenario) ---
#   # This is the **critical new step** that links your clustering to occuN.
#   # It creates the M x N site-to-cell mapping matrix (w).
#   # M = number of sites from clustering
#   # N = number of cells in the raster
#   # (w[i, j] = 1 if cell j is part of site i, 0 otherwise)
#   print(paste("  [Scenario", i, "] 2. Creating W matrix..."))
#   w_matrix_info <- create_w_matrix(clustered_checklists, OR_normalized, cellCovs_df)
#   w_matrix <- w_matrix_info$w              # The M x N sparse matrix
#   M_sites <- nrow(w_matrix)
  
#   # Loop for each of the 100 simulations (stochastic replicates)
#   for (sim_num in 1:N_SIMULATIONS) {
    
#     if (sim_num %% 10 == 0) {
#       print(paste("  [Scenario", i, "] Simulation", sim_num, "of", N_SIMULATIONS, "..."))
#     }
    
#     # --- 5.4. SIMULATE OCCUPANCY & OBSERVATIONS ---
#     # This function uses the `occuN_demo.R` logic:
#     # 1. Calculates cell-level lambda_j (abundance) from cellCovs_df & occ_params
#     # 2. Aggregates to site-level lambda_tilde_i = w %*% lambda_j
#     # 3. Calculates site-level psi_i = 1 - exp(-lambda_tilde_i)
#     # 4. Simulates true site state Z_true (M x 1 vector)
#     # 5. Simulates 'species_observed' for each *checklist* based on Z_true
#     #    and the checklist's obsCovs & det_params.
#     sim_results <- generate_simulated_observations(
#       clustered_checklists = clustered_checklists,
#       cellCovs_df = cellCovs_df,
#       w_matrix = w_matrix,
#       occ_params = occ_params,
#       det_params = det_params,
#       occ_cov_names = occ_cov_names
#     )
#     simulated_checklists_with_obs <- sim_results$checklists
#     Z_true <- sim_results$Z_true # The (M x 1) vector of true occupancy
    
    
#     # --- 5.5. CREATE FULL UNMARKED FRAME (Done ONCE per sim) ---
#     # This function takes the simulated checklist data and:
#     # 1. Pivots it into the M x J observation matrix 'y'
#     # 2. Formats the 'obsCovs' into a list of M x J matrices
#     # 3. Creates the 'siteCovs' data frame (if any)
#     # 4. Bundles everything into an unmarkedFrameOccuN
#     umf_all <- create_unmarked_frame(
#       checklists = simulated_checklists_with_obs,
#       cellCovs = cellCovs_df, # Full N x K cell covariate data frame
#       w = w_matrix,             # Full M x N mapping matrix
#       det_cov_names = det_cov_names,
#       site_ids = w_matrix_info$site_ids # The list of site IDs (1:M)
#     )
    
    
#     # Loop for each of the 25 test-set repeats
#     for (test_rep in 1:N_TEST_REPEATS) {
      
#       # Use a tryCatch to ensure one failed fit doesn't stop the whole script
#       tryCatch({
        
#         # --- 5.6. SPLIT TRAIN/TEST (by site) ---
#         site_indices <- 1:M_sites
#         train_sites <- sample(site_indices, size = floor(0.75 * M_sites))
#         test_sites <- setdiff(site_indices, train_sites)
        
#         # Subset the full unmarkedFrame to create train/test sets
#         umf_train <- umf_all[train_sites, ]
        
#         # Create a test data packet for the evaluation function
#         test_data <- list(
#           umf_test = umf_all[test_sites, ],
#           Z_true_test = Z_true[test_sites]
#         )
        
#         # --- 5.7. FIT occuN MODEL ---
#         # This function fits the occuN model.
#         # It should *internally* run the N_REPEATS=30 random starts
#         # and return only the best-converged model.
#         # (Logic from occuN_demo.R, lines 206-247)
#         model_fit <- fit_occuN(
#           formula = full_formula,
#           data = umf_train,
#           n_starts = N_REPEATS # Use N_REPEATS for random starts
#         )
        
#         # Handle fit failure
#         if (is.null(model_fit)) {
#           stop("Model fit failed for all N_REPEATS starts.")
#         }
        
#         # --- 5.8. EVALUATE MODEL ---
#         # This function uses the fitted model to predict on 'test_data$umf_test',
#         # compares predictions to 'test_data$Z_true_test',
#         # and returns a list of metrics (e.g., auc, rmse).
#         metrics <- calculate_metrics(model_fit, test_data)
        
#         # --- 5.9. STORE SUMMARIZED RESULTS ---
#         result_row <- data.frame(
#           scenario_id = i,
#           clustering_method = clustering_method,
#           simulation_num = sim_num,
#           test_repeat_num = test_rep,
#           auc = metrics$auc,
#           rmse = metrics$rmse,
#           brier = metrics$brier
#           # ... other metrics ...
#         )
        
#         # Add parameter info
#         # This adds columns like "occ_intercept", "elevation", etc.
#         result_row <- cbind(result_row, scenario_row[, c("occ_intercept", occ_cov_names, "det_intercept", det_cov_names)])
        
#         all_results[[length(all_results) + 1]] <- result_row
        
#       }, error = function(e) {
#         # --- 5.10. LOG ERRORS ---
#         print(paste0("  ERROR in scenario=", i, " sim=", sim_num, " test_rep=", test_rep, ": ", e$message))
#         # Store NA results so we know this run failed
#         result_row <- data.frame(
#           scenario_id = i,
#           clustering_method = clustering_method,
#           simulation_num = sim_num,
#           test_repeat_num = test_rep,
#           auc = NA,
#           rmse = NA,
#           brier = NA
#         )
#         result_row <- cbind(result_row, scenario_row[, c("occ_intercept", occ_cov_names, "det_intercept", det_cov_names)])
#         all_results[[length(all_results) + 1]] <- result_row
#       }) # End tryCatch
      
#     } # End test repeat loop
    
#   } # End simulation loop
  
#   # Optional: Save intermediate results after each scenario
#   # final_results_df_temp <- do.call(rbind, all_results)
#   # write_csv(final_results_df_temp, paste0("simulation_experiments/results/temp_results_scenario_", i, ".csv"))
  
# } # End scenario loop

# # -----------------------------------------------
# # 6. SAVE FINAL RESULTS
# # -----------------------------------------------
# print("--- 6. Saving Final Results ---")

# # Combine the list of all results into one big data frame
# final_results_df <- do.call(rbind, all_results)

# # Create results directory if it doesn't exist
# results_dir <- "simulation_experiments/results"
# if (!dir.exists(results_dir)) {
#   dir.create(results_dir, recursive = TRUE)
# }

# # Save the *one* final summary file
# write_csv(final_results_df, file.path(results_dir, "simulation_summary.csv"))

# print(paste("Simulation complete. Final results saved to", file.path(results_dir, "simulation_summary.csv")))