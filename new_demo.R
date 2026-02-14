# new_demo.R
# Purpose: Isolate why occuN is returning NULL/Errors
# (CORRECTED: Fixes Coordinate Projection Error)

source("R/utils.R")
source("R/data_helpers.R")
source("R/clustering_helpers.R")
source("R/model_helpers.R")
source("R/analysis_helpers.R")

library(dplyr)
library(terra)
library(unmarked)
library(sf)

# 1. SETUP DATA (AMCR)
species_name <- "AMCR"
cat(sprintf(">>> DEBUGGING SPECIES: %s\n", species_name))

# Load Config
res_m <- 100 
buffer_m <- 200
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
state_cov_names <- c("elevation", "TCB", "TCG", "TCW", "TCA")
obs_cov_names <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")

# Load & Process Raster
cat("--- Loading Raster ---\n")
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster_raw) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster_raw) <- state_cov_names
cov_tif_albers <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = res_m)
cov_tif_albers <- standardize_state_covs(cov_tif_albers)$raster

# Load Checklists
cat("--- Loading Checklists ---\n")
train_obs <- read.csv(file.path("checklist_data", "species", species_name, paste0(species_name, "_zf_filtered_region_2017.csv")))
train_obs <- train_obs[!is.na(train_obs$duration_minutes),]

# Filter NAs in coordinates immediately
train_obs <- train_obs[!is.na(train_obs$latitude) & !is.na(train_obs$longitude), ]

# 2. RUN CLUSTERING (Kappa=30)
# We calculate clustering on PROJECTED points, but we KEEP original lat/lon in the dataframe
cat("--- Running Clustering (Kappa=30) ---\n")

# Project for K-Means ONLY
train_sf <- st_as_sf(train_obs, coords=c("longitude", "latitude"), crs=4326) %>% st_transform(albers_crs_str)
coords_proj <- st_coordinates(train_sf)

# Prepare dataframe for clustering
train_locs <- train_obs %>% 
  dplyr::select(locality_id, latitude, longitude) %>%
  distinct(locality_id, .keep_all = TRUE)

# We need to match the projected coords to the unique locations
# (Simplest way: project the unique locations)
locs_sf <- st_as_sf(train_locs, coords=c("longitude", "latitude"), crs=4326) %>% st_transform(albers_crs_str)
locs_coords <- st_coordinates(locs_sf)

set.seed(123)
k_val <- round(nrow(train_locs) * 0.30)
km <- kmeans(locs_coords, centers=k_val)
train_locs$site <- km$cluster

# Join site IDs back to main data
# IMPORTANT: clust_df must retain WGS84 latitude/longitude columns for create_site_geometries
clust_df <- train_obs %>%
  dplyr::select(-any_of("site")) %>%
  inner_join(train_locs[, c("locality_id", "site")], by = "locality_id")

# 3. GENERATE W MATRIX
cat("--- Generating W Matrix ---\n")
# create_site_geometries expects WGS84 lat/long, which clust_df now correctly has
current_geoms <- create_site_geometries(clust_df, cov_tif_albers, buffer_m, "Debug_Method")
current_geoms <- st_make_valid(current_geoms) 

w_matrix <- generate_overlap_matrix(current_geoms, cov_tif_albers)
cat(sprintf("    W Matrix Dim: %d x %d\n", nrow(w_matrix), ncol(w_matrix)))
cat(sprintf("    W Matrix Sum: %.2f (Should be > 0)\n", sum(w_matrix)))

# 4. PREPARE UMF
cat("--- Preparing UMF ---\n")
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))[, state_cov_names, drop = FALSE]
full_raster_covs[is.na(full_raster_covs)] <- 0

umf <- prepare_occuN_data(train_obs, clust_df, w_matrix, obs_cov_names, full_raster_covs)
cat("    UMF Summary:\n")
# Print dimensions to check validity
print(dim(umf@y))

# 5. ATTEMPT FIT (With Verbose Error)
cat("\n--- ATTEMPTING FIT ---\n")
obs_formula <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
state_formula <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))

cat("Formula:\n")
print(paste(paste(deparse(obs_formula), collapse=""), paste(deparse(state_formula), collapse="")))

# Run DIRECTLY (no tryCatch wrapper) to see the error
fm <- unmarked::occuN(
  formula = as.formula(paste(paste(deparse(obs_formula), collapse=""), paste(deparse(state_formula), collapse=""))),
  data = umf,
  starts = rep(0, 12), # Try zero starts first
  se = FALSE,
  method = "nlminb"
)

print(fm)