# -----------------------------------------------------------------
# Simulation for occuN model
#
# -----------------------------------------------------------------

############
# 1. Installation & Libraries
############

# ... (your installation code) ...

library(dplyr)
library(foreach)
library(doParallel) # <-- ADDED
library(terra)      # <-- ADDED (for tmpFiles)
library(unmarked)   # <-- ADDED (for occu)
library(rje)        # <-- ADDED (for expit)
library(PRROC)      # <-- ADDED (for AUC/AUPRC)
library(dggridR)    # <-- ADDED (for subsampling)

# Source all your helper functions
source(file.path("R", "utils.R"))
source(file.path("R", "simulation_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))

set.seed(123) # For reproducibility

# ... (your comparison_method_list, optimizer, etc.) ...
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


############
# 1.5. SET UP PARALLEL BACKEND
############
n_cores <- detectCores() - 1 # Use all but one core
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("--- Registered parallel backend with %d cores ---\n", n_cores))

# This is an alternative to using .packages and .export,
# it ensures all workers have the necessary libraries and functions
clusterEvalQ(cl, {
  library(dplyr)
  library(terra)
  library(unmarked)
  library(rje)
  library(PRROC)
  library(dggridR)
  
  # Source all helpers on each worker
  source(file.path("R", "utils.R"))
  source(file.path("R", "simulation_helpers.R"))
  source(file.path("R", "clustering_helpers.R"))
  source(file.path("R", "model_helpers.R"))
  source(file.path("R", "analysis_helpers.R"))
})


###########
# 2. LOAD CONFIGURATION
###########
sim_params <- read.delim(file.path("config","simulation_parameters.csv"), sep = ",", header = T)
sim_clusterings <- read.delim(file.path("config","simulation_clusterings.csv"), sep = ",", header = T)

n_simulations <- 3
n_fit_repeats <- 3 # Note: Your old code uses 'runs' (25). This is for model fitting.
n_test_repeats <- 3  # This is for the test/prediction repeats.

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
    truth_df = NULL
)
cat(sprintf("--- Pre-computing complete. Found %d total clusterings. ---\n", length(all_clusterings)))


cat("--- Pre-computing site geometries for ALL clustering methods... ---\n")
all_site_geometries <- list()
for (method_name in all_method_names) {
    # ... (your geometry creation code) ...
    cat(paste("    - Generating geometries for:", method_name, "\n"))
    
    cluster_data <- all_clusterings[[method_name]]
    if (is.list(cluster_data) && "result_df" %in% names(cluster_data)) {
      cluster_data <- cluster_data$result_df
    }
    
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

# ... (your code to calculate and write clustering_summary_df) ...
clustering_summary_df <- summarize_clusterings(
  all_clusterings = all_clusterings,
  all_site_geometries = all_site_geometries,
  units = "m" 
)
output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(clustering_summary_df, file.path(output_dir, "clustering_descriptive_stats.csv"), row.names = FALSE)
cat(sprintf("--- Summary metrics saved to %s/clustering_descriptive_stats.csv ---\n", output_dir))


# These lists will store results from all simulations
all_dataset_stats <- list()
all_results <- list()


for (cluster_idx in seq_len(nrow(sim_clusterings))) {

  current_clustering_method <- sim_clusterings$method[cluster_idx]
  print(paste("STARTING REFERENCE CLUSTERING:", current_clustering_method))
  
  current_reference_dataframe <- all_clusterings[[current_clustering_method]]
  current_site_geometries <- all_site_geometries[[current_clustering_method]]
  
  cat(paste("    - Pre-calculating Albers projection for:", current_clustering_method, "\n"))
  albers_crs_str <- sf::st_crs(current_site_geometries)$wkt
  cov_tif_albers <- terra::project(state_cov_raster, albers_crs_str, method="bilinear", res = 30)
  area_j_raster <- terra::cellSize(cov_tif_albers, unit="m")
  cat(paste("    - Pre-calculation complete.\n"))

  # Loop over reference parameters (iterating by row index)
  for (param_idx in seq_len(nrow(sim_params))) {

    current_parameter_set <- sim_params[param_idx, ]
    print(paste0("STARTING PARAMETER SET ", param_idx, ":"))
    print(current_parameter_set)
  
    # === 3. REPLACED LOOP ===
    # This is the parallel loop. It replaces `for (sim_num in 1:n_simulations)`
    # We export the key data objects needed by the workers.
    # Functions are available because we used clusterEvalQ.
    sim_results_list <- foreach(
        sim_num = 1:n_simulations,
        .export = c("current_reference_dataframe", "current_site_geometries", 
                    "current_parameter_set", "state_cov_names", "obs_cov_names",
                    "norm_list", "cov_tif_albers", "area_j_raster", "state_cov_raster",
                    "n_test_repeats", "all_clusterings", "selected_optimizer",
                    "current_clustering_method", "param_idx")
    ) %dopar% {
      
      print(paste("--- Simulation run", sim_num, "of", n_simulations, "---"))

      # === 1. SIMULATE DATA (inside worker) ===
      train_data <- simulate_train_data(
        reference_clustering_df = current_reference_dataframe, 
        site_geoms_sf = current_site_geometries,
        parameter_set_row = current_parameter_set, 
        state_cov_names = state_cov_names, 
        obs_cov_names = obs_cov_names,
        norm_list = norm_list,
        cov_tif_albers = cov_tif_albers,
        area_j_raster = area_j_raster
      )
            
      test_data_full <- simulate_test_data(
          norm_list = norm_list, 
          parameter_set_row = current_parameter_set, 
          state_cov_names = state_cov_names, 
          obs_cov_names = obs_cov_names,
          cov_tif = state_cov_raster
      )

      # === 1.5. CALCULATE DATASET STATS (inside worker) ===
      current_dataset_stats <- summarize_datasets(train_data, test_data_full)
      current_dataset_stats$cluster_method <- current_clustering_method
      current_dataset_stats$param_set <- param_idx
      current_dataset_stats$sim_num <- sim_num
      current_dataset_stats <- current_dataset_stats %>%
        dplyr::relocate(cluster_method, param_set, sim_num)
      
      # This list will hold all metrics for THIS sim_num
      sim_run_metrics <- list() 
      model_list <- list()

      # === 3. RUN EXPERIMENTS FOR EACH TEST REPEAT (inside worker) ===
      for (repeat_num in 1:n_test_repeats) {

        # # Spatially subsample the test data
        # set.seed(repeat_num) # Reproducible subsample
        # hexagons <- dggridR::dgconstruct(spacing = 5, topology = "HEXAGON")
        # test_data_full$cell <- dggridR::dgGEO_to_SEQNUM(hexagons, test_data_full$latitude, test_data_full$longitude)$seqnum
        
        # test.det.df <- test_data_full[test_data_full$species_observed == T,]
        # test.undet.df <- test_data_full[test_data_full$species_observed == F,]
        
        # cell_names <- names(table(test.det.df$cell))
        # det.valid.df <- spatial_subsample(test.det.df, cell_names)
        
        # undet.cell_names <- names(table(test.undet.df$cell))
        # undet.valid.df <- spatial_subsample(test.undet.df, undet.cell_names)
        
        # if(nrow(undet.valid.df) > nrow(det.valid.df)){
        #     idx <- sample(seq(1:nrow(undet.valid.df)), nrow(det.valid.df))
        #     undet.valid.df <- undet.valid.df[idx,]
        # } else if (nrow(undet.valid.df) < nrow(det.valid.df)){
        #     idx <- sample(seq(1:nrow(det.valid.df)), nrow(undet.valid.df))
        #     det.valid.df <- det.valid.df[idx,]
        # }
        # test.df <- rbind(det.valid.df, undet.valid.df)

        # # === 4. LOOP OVER EACH CLUSTERING METHOD (inside worker) ===
        for(method_name in names(all_clusterings)){
        
        #   current_clustering_df <- all_clusterings[[method_name]]
        #   set.seed(1) # Consistent seed for modeling
          
        #   # Handle BayesOpt special case
        #   if (is.list(current_clustering_df) && "result_df" %in% names(current_clustering_df)) {
        #       current_clustering_df <- current_clustering_df$result_df
        #   }

        #   # --- Fit models only on the first test repeat ---
        #   if(repeat_num == 1){
            
        #     # (We skip writing files here, as that's very slow in parallel)
            
        #     # Calculate clustering stats (ARI, etc.) against the ground truth
        #     cl_stats <- list(ari=NA, ami=NA , nid=NA)
        #     if (!(method_name %in% c("1to10", "2to10", "2to10-sameObs", "1-per-UL"))){
        #         cl_stats <- calcClusteringStats(current_clustering_df, train_data)
        #     }

        #     # Fit the occupancy model
        #     test.formula <- calcOccModel(current_clustering_df, state_cov_names, obs_cov_names) 
            
        #     occ_par_list <- test.formula@estimates@estimates$state@estimates 
        #     det_par_list <- test.formula@estimates@estimates$det@estimates

        #     model_list[[method_name]] <- list(
        #       occu_parameters = occ_par_list,
        #       det_parameters = det_par_list,
        #       cl_stats = cl_stats,
        #       neg_log_like = test.formula@negLogLike
        #     )
        #   } # --- End of (repeat_num == 1) block ---

        #   # Get the model parameters
        #   model_fit <- model_list[[method_name]]
        #   occ_par_list <- model_fit$occu_parameters
        #   det_par_list <- model_fit$det_parameters

        #   # Calculate predictions on the *current* test.df subsample
        #   test.df$occupied_prob_est <- calculate_weighted_sum(occ_par_list, test.df)
        #   test.df$occupied_prob_est <- rje::expit(test.df$occupied_prob_est)
          
        #   test.df$det_prob_est <- calculate_weighted_sum(det_par_list, test.df)
        #   test.df$det_prob_est <- rje::expit(test.df$det_prob_est)
          
        #   pred_observ <- unlist(test.df$occupied_prob_est * test.df$det_prob_est)

        #   dets <- pred_observ[test.df$species_observed == TRUE]
        #   nondets <- pred_observ[test.df$species_observed == FALSE]

        #   # ROC Curve    
        #   roc <- roc.curve(scores.class0 = dets, scores.class1 = nondets, curve = F)
        #   # PR Curve
        #   pr <- pr.curve(scores.class0 = dets, scores.class1 = nondets, curve = F)
          
        #   occ.mape.i <- mean(abs((test.df$occupied_prob - test.df$occupied_prob_est)/test.df$occupied_prob)) * 100
        #   det.mape.i <- mean(abs((test.df$det_prob - test.df$det_prob_est)/test.df$det_prob)) * 100

          # metrics_df <- data.frame(
          #   cluster_method_ref = current_clustering_method,
          #   param_set_ref = param_idx,
          #   sim_num = sim_num,
          #   test_repeat_num = repeat_num,
          #   method_comp = method_name,
          #   ari = model_fit$cl_stats$ari,
          #   ami = model_fit$cl_stats$ami,
          #   nid = model_fit$cl_stats$nid,
          #   neg_log_like = model_fit$neg_log_like,
          #   occ_prob_mape = occ.mape.i,
          #   det_prob_mape = det.mape.i,
          #   auc = roc$auc,
          #   auprc = pr$auc.integral
          # )

          metrics_df <- data.frame(
            cluster_method_ref = 0
          )

          
          
          # Append to the list for *this worker's* simulation
          sim_run_metrics[[length(sim_run_metrics) + 1]] <- metrics_df

        } # End loop over comparison methods
      } # End loop over test repeats (repeat_num)

      cat(paste("--- Cleaning up temp files for sim", sim_num, "---\n"))
      terra::tmpFiles(remove = TRUE) 
      gc()
      
      # === 5. RETURN RESULTS (from worker) ===
      # This returns the results for a single sim_num
      return(list(
        stats = current_dataset_stats,
        metrics = dplyr::bind_rows(sim_run_metrics) # Combine all metrics for this sim
      ))

    } # --- END of foreach loop ---

    # === 6. PROCESS RESULTS (on main thread) ===
    # Now, loop through the list of results returned by foreach
    # and append them to the global lists
    cat(sprintf("--- Aggregating results for Param Set %d ---\n", param_idx))
    for (result in sim_results_list) {
      all_dataset_stats[[length(all_dataset_stats) + 1]] <- result$stats
      all_results[[length(all_results) + 1]] <- result$metrics
    }

  } # End parameter loop (param_idx)
} # End clustering loop (cluster_idx)


############
# 7. CLEANUP AND SAVE
############
stopCluster(cl)
cat("--- Parallel cluster stopped ---\n")

print("Done")

# After all loops are done
final_dataset_stats_df <- dplyr::bind_rows(all_dataset_stats)
write.csv(final_dataset_stats_df, 
          file.path(output_dir, "dataset_descr_stats.csv"), 
          row.names = FALSE)
cat(sprintf("--- Dataset descriptive stats saved to %s/dataset_descr_stats.csv ---\n", output_dir))

# Bind all simulation metrics
final_results_df <- dplyr::bind_rows(all_results)
write.csv(final_results_df, file.path(output_dir, "simulation_summary.csv"), row.names = FALSE)
cat(sprintf("--- All simulation results saved to %s/simulation_summary.csv ---\n", output_dir))