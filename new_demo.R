# debug_occun.R
# Purpose: Isolate why occuN is returning NULL/Errors

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
# Project to Albers
train_sf <- st_as_sf(train_obs, coords=c("longitude", "latitude"), crs=4326) %>% st_transform(albers_crs_str)
train_obs$latitude <- st_coordinates(train_sf)[,2]
train_obs$longitude <- st_coordinates(train_sf)[,1]

# 2. RUN CLUSTERING (Mock Kappa = 30)
cat("--- Running Clustering (Kappa=30) ---\n")
train_locs <- train_obs %>% distinct(locality_id, latitude, longitude, .keep_all = TRUE)
# Simple K-means for debug speed
set.seed(123)
k_val <- round(nrow(train_locs) * 0.30)
km <- kmeans(train_locs[, c("latitude", "longitude")], centers=k_val)
train_locs$site <- km$cluster

clust_df <- train_obs %>%
  select(-any_of("site")) %>%
  inner_join(train_locs[, c("locality_id", "site")], by = "locality_id")

# 3. GENERATE W MATRIX
cat("--- Generating W Matrix ---\n")
current_geoms <- create_site_geometries(clust_df, cov_tif_albers, buffer_m, "Debug_Method")
# Ensure valid geometries
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
print(summary(umf))

# 5. ATTEMPT FIT (With Verbose Error)
cat("\n--- ATTEMPTING FIT ---\n")
obs_formula <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
state_formula <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))

# Run DIRECTLY (no tryCatch wrapper) to see the error
fm <- unmarked::occuN(
  formula = as.formula(paste(paste(deparse(obs_formula), collapse=""), paste(deparse(state_formula), collapse=""))),
  data = umf,
  starts = rep(0, 12), # Try zero starts first
  se = FALSE,
  method = "nlminb"
)

print(fm)