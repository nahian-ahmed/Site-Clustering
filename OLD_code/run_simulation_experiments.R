##############################################
# Main file for running simulation experiments

# June 27, 2025
##############################################


library(here)

# set working directory to current/source directory
setwd(here::here())

# Import required files
source("helper/BayesOptHelper.R")
source("helper/clustGeoHelper.R")
source("helper/DBSCHelper.R")
source("helper/helpers.R")
source("helper/kmsq.R")
source("helper/mainHelper.R")




#######
# Main function for running site clustering experiments
#######
run_site_clustering_experiments = function(
    species_df,
    method_names,
    occu_model_types,
    runs,
    occ_covs,
    det_covs,
    occ_covs_all,
    det_covs_all,
    simulate=TRUE,
    preprocess=TRUE,
    experiments=TRUE,
    summarize=TRUE,
    plot_results=TRUE) {
    
    species_df_colnames <- c("species", "ref_method", "occ_intercept", occ_covs, "det_intercept", det_covs)
    species_df <- species_df[, species_df_colnames]

    assign("species_df", species_df, envir = .GlobalEnv)
    assign("method_names", method_names, envir = .GlobalEnv)
    assign("species_names", as.list(species_df$species), envir = .GlobalEnv)
    assign("runs", runs, envir = .GlobalEnv)
    assign("occ_covs", occ_covs, envir = .GlobalEnv)
    assign("det_covs", det_covs, envir = .GlobalEnv)
    assign("occ_covs_all", occ_covs_all, envir = .GlobalEnv)
    assign("det_covs_all", det_covs_all, envir = .GlobalEnv)


    if(simulate){
        source("simulation_experiments/simulate_data.R")
    }
    assign("species_df_formatted", read.delim("results/simulation_experiments/summarized/simulation_species_formatted.csv", sep = ",", header = T), envir = .GlobalEnv)    

    if (preprocess){
        source("simulation_experiments/preprocess_data.R")
    }


    if (experiments){
        
        for (species_name in species_names){

            assign("spec_name", species_name, envir = .GlobalEnv)

            # Run experiments on specific simulated species data
            set.seed(1) # Ensures that order of species in species_names does not affect reproducibility
            source("simulation_experiments/run_experiments.R")

        }
    }

    if (summarize){
        # Summarize results
        source("simulation_experiments/summarize_results.R")
    }

    if (plot_results){
        # Plot results
        source("simulation_experiments/plot_results.R")
    }
}




#######
# Specify simulation species to run experiments on
# Please see simulation_species.csv for more information
#######
species_df <- read.delim("simulation_species.csv", sep = ",", header = T)


#######
# Specify site clustering methods
# Please see helper/mainHelper.R for method strings
#######
method_names <- c(
    "reference_clustering",
    "one_to_10",
    "two_to_10",
    "two_to_10_sameObs",
    "kmSq-1000",
    "lat_long",
    "rounded-4",
    "svs",
    "one_UL",
    "DBSC",
    "BayesOptClustGeo"
)



#######
# Specify number of experiment runs/repeats
#######
runs <- 25

#######
# Specify occupancy and detection covariates
#######
# elevation and time_observations_started
occ_covs <- c("elevation", "TCB", "TCG", "TCW", "TCA")
det_covs <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")


#######
# Run site clustering experiments
#######
run_site_clustering_experiments(
    species_df = species_df, 
    method_names = method_names,
    runs = runs,
    occ_covs = occ_covs,
    det_covs = det_covs,
    occ_covs_all = occ_covs,
    det_covs_all = det_covs,
    simulate = FALSE, 
    preprocess = FALSE, 
    experiments = FALSE, 
    summarize = FALSE, 
    plot_results = TRUE
)
