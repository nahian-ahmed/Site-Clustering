# -----------------------------------------------------------------
# Species Point-Based Experiments
# Comparison: Buffered vs Unbuffered
# -----------------------------------------------------------------

###
# 1. SETUP
###

install_now = FALSE
if (install_now){
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  if (!requireNamespace("devtools", quietly = FALSE)) install.packages("devtools")
  suppressMessages(devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE))
}

library(unmarked)
library(dplyr)
library(tidyr)
library(terra)
library(PRROC)
library(sf)

# Source Helpers
source(file.path("R", "utils.R"))
source(file.path("R", "data_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))

# --- CONFIGURATION ---
species_names <- c(
  "AMCR", "AMRO", "BAEA", "BKHGRO", "BRCR", "BUTI", "CASC", "CHBCHI", 
  "COHA", "HAFL", "HAWO", "HEWA", "MAWA", "MOQU", "NOFL", "NOOW", 
  "OLFL", "PAFL", "PAWR", "PIWO", "REHA", "SOSP", "SPTO", "SWTH", 
  "WAVI", "WEPE", "WETA", "WIWA", "WRENTI", "YEBCHA", "YEWA"
)
method_names <- c("lat-long", "1to10", "2to10")


# Covariates
state_cov_names <- c("elevation", "TCB", "TCG", "TCW", "TCA")
obs_cov_names <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")

# Optimization & Simulation Settings
selected_optimizer <- "nlminb"
n_fit_repeats <- 50
n_test_repeats <- 25

res_m <- 100 
buffers_m  <- c(100, 200, 500)
test_buffer_m <- 200
hex_m <- 100



# Output
output_dir <- file.path("species_experiments", "point_based_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)








# Save
final_res <- bind_rows(results_list)
write.csv(final_res, file.path(output_dir, "point_based_results.csv"), row.names=FALSE)
cat("Done.\n")