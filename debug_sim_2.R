# -----------------------------------------------------------------
# debug_sim_2.R
# Diagnostic script to debug Parameter Recovery failure
# -----------------------------------------------------------------

library(unmarked)
library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(Matrix)
library(rje)

# Source existing helpers
source(file.path("R", "utils.R"))
source(file.path("R", "simulation_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))
source(file.path("R", "plotting_helpers.R"))

set.seed(123)

# 1. SETUP SINGLE VARIANT (V1)
cat("\n=== 1. SETUP DEBUG RUN (V1_Uniform_Simple) ===\n")
sim_params <- read.csv(file.path("config", "simulation_parameters.csv"))
# Pick the first parameter set
curr_params <- sim_params[1, ]
cat("True State Intercept from Config:", curr_params$state_intercept, "\n")

# Load Raster
state_cov_names_all <- names(sim_params)[2:6] 
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster_raw) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster_raw) <- state_cov_names_all
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, 
                                     "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", 
                                     method="bilinear", res = 100)
standardization_results <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- standardization_results$raster
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))
full_raster_covs[is.na(full_raster_covs)] <- 0

# 2. GENERATE DATA (Small Batch)
cat("\n=== 2. GENERATING DATA (N=1000) ===\n")
# Create dummy boundary for this test
boundary_vect <- terra::vect(file.path("state_covariate_raster", "boundary", "boundary.shp"))
boundary_vect <- terra::project(boundary_vect, terra::crs(cov_tif_albers))

# Generate Uniform Data
masked <- terra::mask(cov_tif_albers[[1]], boundary_vect)
valid_cells <- terra::cells(masked)
sampled_idx <- sample(valid_cells, 1000, replace = TRUE)
coords_proj <- terra::xyFromCell(cov_tif_albers[[1]], sampled_idx)
df_proj <- data.frame(x = coords_proj[,1], y = coords_proj[,2])
v_proj <- terra::vect(df_proj, geom=c("x", "y"), crs = terra::crs(cov_tif_albers))
v_geo <- terra::project(v_proj, "+proj=longlat +datum=WGS84")
coords_geo <- terra::crds(v_geo)

base_df <- data.frame(
  checklist_id = paste0("debug_", 1:1000),
  locality_id = paste0("loc_debug_", 1:1000),
  latitude = coords_geo[,2],
  longitude = coords_geo[,1],
  observation_date = "2017-06-01",
  formatted_date = "2017-06-01",
  duration_minutes = rnorm(1000)
)
env_df <- extract_state_covs(base_df, cov_tif_albers)
base_df <- dplyr::inner_join(base_df, env_df, by = "checklist_id")

# 3. CLUSTERING (2-kmSq)
cat("\n=== 3. RUNNING CLUSTERING (2-kmSq) ===\n")
cluster_res <- run_clustering_method("2-kmSq", base_df, names(state_cov_raster_raw))
cluster_df <- cluster_res$data

# Geometries
geoms <- create_site_geometries(cluster_df, cov_tif_albers, 200, "2-kmSq")
split_res <- disjoint_site_geometries(geoms, cluster_df)
final_geoms <- split_res$geoms
final_cluster_df <- split_res$data

# W Matrix
w_matrix <- generate_overlap_matrix(final_geoms, cov_tif_albers)

cat("W Matrix Stats:\n")
cat("  Dims:", dim(w_matrix), "\n")
cat("  Sum of entire W:", sum(w_matrix), "km2\n")
cat("  Mean Area per Site:", mean(rowSums(w_matrix)), "km2\n")
if(mean(rowSums(w_matrix)) < 0.01) warning("!!! WARNING: W matrix entries are extremely small. Unit mismatch? !!!")

# 4. SIMULATION
cat("\n=== 4. SIMULATING ABUNDANCE & OBSERVATIONS ===\n")
# Calculate Density
state_intercept <- curr_params$state_intercept
log_lambda_j <- cov_tif_albers[[1]] * 0 + state_intercept
# Add elevation if used
log_lambda_j <- log_lambda_j + (cov_tif_albers[["elevation"]] * curr_params$elevation)

cell_density_val <- terra::values(exp(log_lambda_j), mat = FALSE)
cell_density_val[is.na(cell_density_val)] <- 0

cat("Simulating with:\n")
cat("  State Intercept:", state_intercept, "\n")
cat("  Mean Cell Density (raw):", mean(cell_density_val), "\n")

# Expected Abundance
lambda_tilde_i <- as.numeric(w_matrix %*% cell_density_val)
cat("  Mean Site Abundance (lambda_tilde):", mean(lambda_tilde_i), "\n")

# Simulate N
N <- rpois(length(lambda_tilde_i), lambda_tilde_i)
cat("  Mean Simulated N:", mean(N), "\n")
cat("  Occupancy Rate:", mean(N > 0), "\n")

# Simulate Observations
obs_par_list <- list(intercept = curr_params$obs_intercept, duration_minutes = 0) # simplified
train_data <- simulate_train_data(
  final_cluster_df, 
  c("duration_minutes"), 
  obs_par_list, 
  w_matrix, 
  cell_density_val
)

det_count <- sum(train_data$species_observed)
cat("  Total Detections in Data:", det_count, "\n")
if(det_count == 0) stop("!!! STOP: No detections simulated. Params too low? !!!")

# 5. DATA PREP FOR MODEL
cat("\n=== 5. PREPARING occuN DATA ===\n")
umf <- prepare_occuN_data(
  train_data, 
  final_cluster_df, 
  w_matrix, 
  c("duration_minutes"), 
  full_raster_covs
)

y_mat <- umf@y
cat("Y Matrix Stats:\n")
cat("  Dims:", dim(y_mat), "\n")
cat("  Number of sites with at least one detection:", sum(rowSums(y_mat, na.rm=T) > 0), "\n")
cat("  Number of rows that are all NA:", sum(apply(y_mat, 1, function(x) all(is.na(x)))), "\n")

if(sum(rowSums(y_mat, na.rm=T) > 0) == 0) {
  stop("!!! STOP: Y matrix has 0 detections. Join failed or IDs mismatch? !!!")
}

# 6. FITTING
cat("\n=== 6. FITTING MODEL (With Trace) ===\n")
# Start close to truth to check if it holds
truth_starts <- c(
  obs_intercept = curr_params$obs_intercept, 
  obs_dur = 0, 
  state_intercept = curr_params$state_intercept, 
  state_elev = curr_params$elevation
)

cat("Fitting with starts at Truth:\n")
print(truth_starts)

fm <- unmarked::occuN(
  formula = ~duration_minutes ~elevation,
  data = umf,
  starts = truth_starts,
  control = list(trace = 1, REPORT = 1) # Print optimization steps
)

cat("\n=== FINAL ESTIMATES ===\n")
print(coef(fm))

cat("\nDone.\n")