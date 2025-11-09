# -----------------------------------------------
# 1. SETUP
# -----------------------------------------------
# Load libraries

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install forked 'unmarked' package
suppressMessages(
    devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE, quiet = FALSE)
)


library(unmarked) # (Your custom occuN version)
library(tidyverse)
# ... other libraries

# Source all your new helper functions
source("R/utils.R")
source("R/simulation_helpers.R")
source("R/clustering_helpers.R")
source("R/model_helpers.R")

# -----------------------------------------------
# 2. LOAD CONFIGURATION
# -----------------------------------------------
# Load the 81 "species" scenarios
sim_params <- read_csv("config/simulation_params.csv")
sim_clusterings <- read_csv("config/simulation_clusterings.csv")

# Combine them into one 81-row "scenario" data frame
scenarios <- cross_join(sim_params, sim_clusterings)

# Define simulation loops
N_SIMULATIONS <- 100 # "multiple simulations (100)"
N_REPEATS <- 30     # "repeats (30)"
N_TEST_REPEATS <- 25  # "test set repeats (25)"
# Note: You'll need to decide how these loops are nested.
# Below is one possible interpretation.

# -----------------------------------------------
# 3. INITIALIZE RESULTS
# -----------------------------------------------
# Create an empty list to store the final metrics from every run.
# Storing a data frame is fine too, but a list is often more flexible.
all_results <- list()

# -----------------------------------------------
# 4. RUN MASTER LOOP
# -----------------------------------------------
# Loop over each of the 81 scenarios
for (i in 1:nrow(scenarios)) {
  
  scenario_row <- scenarios[i, ]
  
  # Loop for each of the 100 simulations
  for (sim_num in 1:N_SIMULATIONS) {
    
    # === THIS IS THE "IN-ONE-GO" PIPELINE ===
    
    # 1. SIMULATE
    # Generates truth (with variable site areas) + observations
    raw_data <- generate_simulation_data(scenario_row)
    
    # 2. PREPROCESS / CLUSTER
    # Applies the clustering method for this scenario
    clustered_data <- apply_clustering(raw_data, scenario_row$clustering_method)
    
    # 3. RUN EXPERIMENT (with test repeats)
    # This loop handles your 25 test-set repeats (e.g., cross-validation)
    for (test_rep in 1:N_TEST_REPEATS) {
      
      # a. Split into train/test
      split_data <- create_train_test_split(clustered_data)
      
      # b. Fit Model (using your new helper)
      # This is where occuN is finally called!
      model_fit <- fit_occuN(
        formula = ~cov1 ~cov2, 
        data = split_data$train,
        site_areas = split_data$train$areas
      )
      
      # c. Evaluate Model
      metrics <- calculate_metrics(model_fit, split_data$test)
      
      # d. Store SUMMARIZED results (not the model object!)
      result_row <- data.frame(
        scenario_id = i,
        simulation_num = sim_num,
        test_repeat_num = test_rep,
        rmse = metrics$rmse,
        auc = metrics$auc
        # ... other metrics
      )
      
      all_results[[length(all_results) + 1]] <- result_row
      
    } # End test repeat loop
    
    # (The 'N_REPEATS <- 30' loop might fit in here or around the
    # simulation loop, depending on your experimental design)
    
    # === END OF "IN-ONE-GO" PIPELINE ===
    
  } # End simulation loop
  
  # Optional: Print progress
  print(paste("Finished scenario", i, "of", nrow(scenarios)))
  
} # End scenario loop

# -----------------------------------------------
# 5. SAVE FINAL RESULTS
# -----------------------------------------------
# Combine the list of results into one big data frame
final_results_df <- do.call(rbind, all_results)

# Save the *one* final summary file
write_csv(final_results_df, "simulation_experiments/results/simulation_summary.csv")

print("Simulation complete. Final results saved.")