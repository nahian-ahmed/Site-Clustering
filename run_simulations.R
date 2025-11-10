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



###########
# 2. LOAD CONFIGURATION
###########
sim_params <- read.csv(file.path("config","simulation_parameters.csv"))
sim_clusterings <- read.csv(file.path("config","simulation_clusterings.csv"))


n_simulations <- 50 
n_fit_repeats <- 25
n_test_repeats <- 25



state_cov_names <- names(sim_params)[2:6]
obs_cov_names <- names(sim_params)[8:12]


# optimizers = "BFGS", "L-BFGS-B", "CG", "Nelder-Mead", "SANN", "nlminb" 
selected_optimizer <- "nlminb"


# Load the base landscape raster
state_cov_raster <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
names(state_cov_raster) <- state_cov_names

base_train_data <- prepare_train_data(state_cov_raster, state_cov_names, obs_cov_names)
base_train_df <- base_train_data$train_df
norm_list <- base_train_data$norm_list



all_results <- list()


# Loop over reference clusterings (iterating by row index)
for (cluster_idx in seq_len(nrow(sim_clusterings))) {

  current_clustering_method <- sim_clusterings$method[cluster_idx]
  
  print(current_clustering_method)
  
  # Loop over reference parameters (iterating by row index)
  for (param_idx in seq_len(nrow(sim_params))) {

    current_parameter_set <- sim_params[param_idx, ]

    print(current_parameter_set)
  
    # Loop for each of the 100 simulations (stochastic replicates)
    for (sim_num in 1:n_simulations) {

      
      # Simualate train data and test data    

      train_data <- simulate_train_data()
      test_data <- simulate_test_data()

      # Fit model n_fit_times and pick fit with lowest NLL


      # set random seed inside spatial subsample. 


      for (repeat_num in 1:n_test_repeats) {

        
      
      
      }


        
      

    
    } # End simulation loop
    
  } # End parameter loop
} # End clustering loop
# --- MODIFIED LOOP END ---



output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the final summary file
write.csv(all_results, file.path(output_dir, "simulation_summary.csv"), row.names = FALSE)

