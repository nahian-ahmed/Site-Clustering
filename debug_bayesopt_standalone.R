# debug_bayesopt_standalone.R
# Purpose: Replicate the EXACT crash condition inside BayesOpt

# 1. Source EVERYTHING
source("R/utils.R")
source("R/data_helpers.R")
source("R/clustering_helpers.R")
source("R/model_helpers.R")
source("R/analysis_helpers.R")
source("R/clustering/clustgeo.R")

library(dplyr)
library(terra)
library(unmarked)
library(sf)
library(ClustGeo)

# 2. Config
species_name <- "AMCR"
res_m <- 100 
buffer_m <- 200
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
state_cov_names <- c("elevation", "TCB", "TCG", "TCW", "TCA")
obs_cov_names <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")

cat(">>> SETUP COMPLETE. LOADING DATA...\n")

# 3. Load Raster & Standardize (Exactly like species_clusters.R)
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster_raw) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster_raw) <- state_cov_names
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = res_m)

# Standardize Raster
standardization_results <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- standardization_results$raster
state_cov_params <- standardization_results$params

# 4. Prepare Train Data (This performs the Standardization)
cat(">>> PREPARING TRAIN DATA (With Scaling)...\n")
train_data_res <- prepare_train_data(
    state_covs = state_cov_names, 
    obs_covs = obs_cov_names, 
    cov_tif = cov_tif_albers_raw, 
    state_standardization_params = state_cov_params,
    placeholder_spec_name = species_name
)
master_train_df <- train_data_res$train_df

# Check scaling
cat("    [Check] Day of Year Mean (Should be ~0):", mean(master_train_df$day_of_year), "\n")
cat("    [Check] Time Obs Mean (Should be ~0):", mean(master_train_df$time_observations_started), "\n")

# 5. Mimic the BayesOpt Data Filtering
# (Min 2 points, Max 10 points)
train_df_spatial_input <- master_train_df %>%
  group_by(locality_id) %>%
  filter(n() >= 2) %>% 
  slice_sample(n = 10) %>%
  ungroup()

# Add species obs (Already in master_train_df, but mimicking the loop join logic)
spec_train_obs <- read.csv(file.path("checklist_data", "species", species_name, paste0(species_name, "_zf_filtered_region_2017.csv")))
spec_train_obs <- spec_train_obs[, c("checklist_id", "species_observed")]

# Create current_spatial_input
current_spatial_input <- train_df_spatial_input %>%
    dplyr::select(-any_of("species_observed")) %>%
    inner_join(spec_train_obs, by = "checklist_id")

# 6. RUN ONE ROUND OF CLUSTERING (Kappa = 30)
cat(">>> RUNNING CLUSTGEO (Kappa = 30)...\n")

# Prep distances
train_locs <- current_spatial_input %>%
    dplyr::distinct(locality_id, latitude, longitude, .keep_all = TRUE)
  
df_for_dist <- train_locs
df_for_dist$latitude  <- as.numeric(scale(df_for_dist$latitude))
df_for_dist$longitude <- as.numeric(scale(df_for_dist$longitude))

env_dist <- dist(df_for_dist[, state_cov_names])
geo_dist <- dist(df_for_dist[, c("latitude", "longitude")])

tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = 0.5)
K <- max(2, round(nrow(train_locs) * 0.30))
train_locs$site <- cutree(tree, K)

clust_df <- current_spatial_input %>%
    dplyr::select(-any_of("site")) %>%
    dplyr::inner_join(train_locs[, c("locality_id", "site")], by = "locality_id")

# 7. GEOMETRY GENERATION (The suspected crash site)
cat(">>> GENERATING GEOMETRIES & W MATRIX...\n")

current_geoms <- create_site_geometries(clust_df, cov_tif_albers, buffer_m, "Debug_Standalone")
current_geoms <- st_make_valid(current_geoms) # Ensure validity

# THE SPLIT STEP
split_res <- disjoint_site_geometries(current_geoms, clust_df)
current_geoms <- split_res$geoms
clust_df <- split_res$data 

w_matrix <- generate_overlap_matrix(current_geoms, cov_tif_albers)

cat(sprintf("    W Matrix: %d sites x %d cells\n", nrow(w_matrix), ncol(w_matrix)))
cat(sprintf("    W Matrix Sum: %.2f\n", sum(w_matrix)))

# 8. FIT OCCUN
cat(">>> FITTING OCCUN...\n")
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))[, state_cov_names, drop = FALSE]
full_raster_covs[is.na(full_raster_covs)] <- 0

umf <- prepare_occuN_data(current_spatial_input, clust_df, w_matrix, obs_cov_names, full_raster_covs)

obs_formula <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
state_formula <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))

# verbose = TRUE to see optimizer output
fm <- fit_occuN_model(
    umf, state_formula, obs_formula,
    n_reps = 1, stable_reps = 1,
    optimizer = "nlminb"
)

if (is.null(fm)) {
    cat(">>> MODEL FAILED (NULL Returned)\n")
} else {
    cat(">>> MODEL SUCCESS!\n")
    print(fm)
}