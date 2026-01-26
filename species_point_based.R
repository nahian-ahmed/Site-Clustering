# -----------------------------------------------------------------
# Species Point-Based Experiments
# Comparison: Buffered vs Unbuffered
# -----------------------------------------------------------------

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
species_names <- c("AMCR") # Add your full list here
methods       <- c("lat-long", "1to10", "2to10")
buffer_sizes  <- c(100, 200, 500) # For the "Buffered" experiments
test_buffer   <- 200              # Fixed Test Buffer for "Buffered" experiments
resolution    <- 100              # Fixed Raster Resolution

# Output
output_dir <- file.path("species_experiments", "point_based_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 1. LOAD & PREP RASTER (STATIC) ---
cat("--- Loading and Standardizing Raster (100m) ---\n")
albers_crs <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +units=m"

state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs, res = resolution)

# Standardize
std_res <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- std_res$raster
std_params <- std_res$params

# Area Raster (for occuN)
cell_area_km2 <- (resolution/1000)^2
area_j_raster <- cov_tif_albers[[1]] * 0 + cell_area_km2
names(area_j_raster) <- "area"

full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))
# Replace NAs with 0 (assuming mean centered)
full_raster_covs[is.na(full_raster_covs)] <- 0

# --- HELPER: Aggregation ---
# Helper to aggregate covariates over polygons for "Buffered occu"
extract_mean_covs <- function(geoms_sf, raster_stack) {
  vect_geoms <- terra::vect(geoms_sf)
  # Extract mean values for each polygon
  means <- terra::extract(raster_stack, vect_geoms, fun = mean, na.rm = TRUE, ID = FALSE)
  return(means)
}

# --- 2. PREPARE TEST DATASETS (Dual Modes) ---
cat("--- Preparing Test Data (Dual Modes) ---\n")

# A. Base Data (Point Extraction / Unbuffered)
# prepare_test_data extracts at exact points by default
test_df_unbuffered <- prepare_test_data(
    state_covs = names(cov_tif_albers),
    obs_covs = c("day_of_year", "time_observations_started", "duration_minutes"),
    cov_tif = cov_tif_albers_raw, # Use RAW for extraction, it gets standardized inside
    standardization_params = std_params,
    placeholder_spec_name = species_names[1]
)

# B. Buffered Spatial Structures (Fixed 200m)
# This creates geometries and W matrix for occuN evaluation
test_structs_buffered <- prepare_test_spatial_structures(
    test_df = test_df_unbuffered, # Starts with df, adds geoms
    albers_crs = albers_crs,
    buffer_m = test_buffer, # FIXED at 200m
    cov_raster_albers = cov_tif_albers,
    area_raster = area_j_raster
)
test_df_buffered <- test_structs_buffered$test_df
test_w_matrix    <- test_structs_buffered$w_matrix

# C. Buffered Covariates for OCCU (Mean over 200m)
# For fair comparison, Buffered OCCU should see the mean habitat of the buffer, not the point.
cat("  - Extracting mean test covariates for Buffered occu...\n")
test_geoms_vect <- terra::vect(voronoi_clipped_buffers(
    sf::st_as_sf(test_df_unbuffered, coords=c("longitude","latitude"), crs=4326) %>% sf::st_transform(albers_crs), 
    buffer_dist = test_buffer
))
test_means <- terra::extract(cov_tif_albers, test_geoms_vect, fun=mean, na.rm=TRUE, ID=FALSE)
test_df_buffered_occu <- cbind(test_df_buffered, test_means) # Add mean cols (will overwrite or duplicate, handled later)


# --- 3. MAIN EXPERIMENT LOOP ---
results_list <- list()

for (sp in species_names) {
  cat(sprintf("\n>>> SPECIES: %s\n", sp))
  
  # Load Obs Data
  train_obs <- read.csv(file.path("checklist_data", "species", sp, paste0(sp, "_zf_filtered_region_2017.csv")))
  test_obs  <- read.csv(file.path("checklist_data", "species", sp, paste0(sp, "_zf_filtered_region_2018.csv")))
  train_obs <- train_obs[!is.na(train_obs$duration_minutes),]

  # Update Test DFs with current species obs
  # 1. Unbuffered
  curr_test_unbuf <- test_df_unbuffered
  curr_test_unbuf$species_observed <- NULL
  curr_test_unbuf <- inner_join(curr_test_unbuf, test_obs[, c("checklist_id", "species_observed")], by="checklist_id")
  
  # 2. Buffered (for occuN - needs W matrix reference)
  curr_test_buf <- test_df_buffered
  curr_test_buf$species_observed <- NULL
  curr_test_buf <- inner_join(curr_test_buf, test_obs[, c("checklist_id", "species_observed")], by="checklist_id")

  # 3. Buffered (for occu - needs Mean Covs)
  curr_test_buf_occu <- test_df_buffered_occu
  curr_test_buf_occu$species_observed <- NULL
  curr_test_buf_occu <- inner_join(curr_test_buf_occu, test_obs[, c("checklist_id", "species_observed")], by="checklist_id")


  # =========================================================
  # EXPERIMENT A: UNBUFFERED (AAAI Style)
  # Buffer = 0 (Points only) | Model = occu ONLY
  # =========================================================
  cat("  [Experiment A] Unbuffered (Point-based, Buffer=0)\n")
  
  # Base Train (Point Extraction)
  base_train_unbuf <- prepare_train_data(
      state_covs = names(cov_tif_albers),
      obs_covs = c("day_of_year", "time_observations_started", "duration_minutes"),
      cov_tif = cov_tif_albers_raw,
      state_standardization_params = std_params,
      placeholder_spec_name = sp
  )$train_df

  for (method in methods) {
      cat(sprintf("    - Method: %s... ", method))
      
      # Filter
      filtered_train <- base_train_unbuf
      if (method == "lat-long") {
        filtered_train$site <- filtered_train$locality_id
      } else if (method == "1to10") {
        filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() <= 10) %>% ungroup()
        filtered_train$site <- filtered_train$locality_id
      } else if (method == "2to10") {
        filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() >= 2, n() <= 10) %>% ungroup()
        filtered_train$site <- filtered_train$locality_id
      }
      
      if(nrow(filtered_train) < 50) { cat("Skip (low N)\n"); next }

      # FIT OCCU (Unbuffered)
      # We use the point-extracted covariates directly found in filtered_train
      tryCatch({
        # Pivot for unmarked (Site x Visit)
        # For Lat-Long/1to10, sites are mostly singletons, but pivot handles repeats if they exist
        occu_df <- filtered_train %>%
          select(site, species_observed, elevation, TCB, TCG, TCW, TCA, 
                 day_of_year, time_observations_started, duration_minutes) %>%
          group_by(site) %>% mutate(visit = row_number()) %>% ungroup()
        
        y_wide <- occu_df %>% select(site, visit, species_observed) %>% 
                  pivot_wider(names_from=visit, values_from=species_observed) %>% select(-site) %>% as.matrix()
        
        # Site Covs (Take first entry per site - assumes point constant)
        site_covs <- occu_df %>% group_by(site) %>% slice(1) %>% 
                     select(elevation, TCB, TCG, TCW, TCA) %>% as.data.frame()
        
        # Obs Covs (Complex pivot, simplified here for robustness)
        obs_covs_list <- list()
        for(v in c("day_of_year", "time_observations_started", "duration_minutes")){
             obs_covs_list[[v]] <- occu_df %>% select(site, visit, all_of(v)) %>% 
                                   pivot_wider(names_from=visit, values_from=all_of(v)) %>% select(-site) %>% as.matrix()
        }
        
        umf <- unmarkedFrameOccu(y = y_wide, siteCovs = site_covs, obsCovs = obs_covs_list)
        
        fm <- occu(~day_of_year + time_observations_started + duration_minutes ~elevation + TCB + TCG + TCW + TCA, umf)
        
        # PREDICT (On Unbuffered Test Data)
        # Using point covariates
        X_state <- model.matrix(~elevation + TCB + TCG + TCW + TCA, curr_test_unbuf)
        X_obs   <- model.matrix(~day_of_year + time_observations_started + duration_minutes, curr_test_unbuf)
        
        beta <- coef(fm, type="state")
        alpha <- coef(fm, type="det")
        
        pred_psi <- plogis(X_state %*% beta)
        pred_det <- plogis(X_obs %*% alpha)
        pred_obs <- pred_psi * pred_det
        
        met <- calculate_classification_metrics(pred_obs, curr_test_unbuf$species_observed)
        
        results_list[[length(results_list)+1]] <- data.frame(
           Species = sp, Buffer = 0, Method = method, Model = "occu",
           AUC = met$auc, AUPRC = met$auprc
        )
      }, error = function(e) cat("Error:", e$message, "\n"))
      
      cat("Done.\n")
  } # End Method Loop (Unbuffered)


  # =========================================================
  # EXPERIMENT B: BUFFERED (occuN Style)
  # Buffer = 100, 200, 500
  # Models = occuN AND occu (using mean covs)
  # =========================================================
  
  for (buf in buffer_sizes) {
    cat(sprintf("  [Experiment B] Buffered (Size: %dm)\n", buf))
    
    # Reload Base Train (Raw - we will extract mean covs manually/via geometries)
    base_train_raw <- prepare_train_data(
        state_covs = names(cov_tif_albers),
        obs_covs = c("day_of_year", "time_observations_started", "duration_minutes"),
        cov_tif = cov_tif_albers_raw,
        state_standardization_params = std_params,
        placeholder_spec_name = sp
    )$train_df
    
    for (method in methods) {
      cat(sprintf("    - Method: %s... ", method))
      
      # Filter
      filtered_train <- base_train_raw
      if (method == "lat-long") {
        filtered_train$site <- filtered_train$locality_id
      } else if (method == "1to10") {
        filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() <= 10) %>% ungroup()
        filtered_train$site <- filtered_train$locality_id
      } else if (method == "2to10") {
        filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() >= 2, n() <= 10) %>% ungroup()
        filtered_train$site <- filtered_train$locality_id
      }
      if(nrow(filtered_train) < 50) { cat("Skip\n"); next }

      # Create Geometries & W Matrix
      train_geoms <- create_site_geometries(filtered_train, cov_tif_albers, buffer_m = buf)
      train_w <- generate_overlap_matrix(train_geoms, cov_tif_albers)
      
      # Extract MEAN covariates for occu
      train_means <- extract_mean_covs(train_geoms, cov_tif_albers)
      train_means$site <- train_geoms$site
      
      # Join Means back to filtered_train for occu
      # (filtered_train currently has point covs, we want to replace or augment)
      # We'll create a specific df for occu fit
      train_df_occu_buf <- filtered_train %>%
        select(-elevation, -TCB, -TCG, -TCW, -TCA) %>% # Drop point covs
        inner_join(train_means, by="site")           # Add mean covs
        

      # --- FIT occuN ---
      tryCatch({
          site_lookup <- filtered_train %>% select(checklist_id, site)
          umf <- prepare_occuN_data(filtered_train, site_lookup, train_w, 
                                    c("day_of_year", "time_observations_started", "duration_minutes"), 
                                    full_raster_covs)
          fm <- fit_occuN_model(umf, ~elevation + TCB + TCG + TCW + TCA, 
                                     ~day_of_year + time_observations_started + duration_minutes, 
                                     n_reps=3, stable_reps=3)
          if(!is.null(fm)) {
             # Predict on Buffered Test Set (Using Area)
             alphas <- coef(fm, 'det'); betas  <- coef(fm, 'state')
             X_state <- model.matrix(~elevation + TCB + TCG + TCW + TCA, curr_test_buf)
             X_obs   <- model.matrix(~day_of_year + time_observations_started + duration_minutes, curr_test_buf)
             
             lambda <- exp(X_state %*% betas)
             psi    <- 1 - exp(-lambda * curr_test_buf$area_j) # <--- AREA DEPENDENT
             p      <- plogis(X_obs %*% alphas)
             
             met <- calculate_classification_metrics(psi * p, curr_test_buf$species_observed)
             results_list[[length(results_list)+1]] <- data.frame(
               Species = sp, Buffer = buf, Method = method, Model = "occuN",
               AUC = met$auc, AUPRC = met$auprc
             )
          }
      }, error = function(e) cat("occuN err: ", e$message, "\n"))


      # --- FIT occu (Buffered) ---
      tryCatch({
          # Use train_df_occu_buf (Mean Covs)
          occu_df <- train_df_occu_buf %>%
            group_by(site) %>% mutate(visit = row_number()) %>% ungroup()
          
          y_wide <- occu_df %>% select(site, visit, species_observed) %>% 
                    pivot_wider(names_from=visit, values_from=species_observed) %>% select(-site) %>% as.matrix()
          
          site_covs <- occu_df %>% group_by(site) %>% slice(1) %>% 
                       select(elevation, TCB, TCG, TCW, TCA) %>% as.data.frame()
          
          obs_covs_list <- list()
          for(v in c("day_of_year", "time_observations_started", "duration_minutes")){
             obs_covs_list[[v]] <- occu_df %>% select(site, visit, all_of(v)) %>% 
                                   pivot_wider(names_from=visit, values_from=all_of(v)) %>% select(-site) %>% as.matrix()
          }
          
          umf <- unmarkedFrameOccu(y = y_wide, siteCovs = site_covs, obsCovs = obs_covs_list)
          fm <- occu(~day_of_year + time_observations_started + duration_minutes ~elevation + TCB + TCG + TCW + TCA, umf)
          
          # Predict on Buffered Test Set (Using Mean Covs)
          # curr_test_buf_occu has the MEAN covariates of the 200m buffer
          X_state <- model.matrix(~elevation + TCB + TCG + TCW + TCA, curr_test_buf_occu)
          X_obs   <- model.matrix(~day_of_year + time_observations_started + duration_minutes, curr_test_buf_occu)
          
          pred_psi <- plogis(X_state %*% coef(fm, "state"))
          pred_det <- plogis(X_obs %*% coef(fm, "det"))
          
          met <- calculate_classification_metrics(pred_psi * pred_det, curr_test_buf_occu$species_observed)
          results_list[[length(results_list)+1]] <- data.frame(
             Species = sp, Buffer = buf, Method = method, Model = "occu",
             AUC = met$auc, AUPRC = met$auprc
          )
      }, error = function(e) cat("occu err: ", e$message, "\n"))
      
      cat("Done.\n")
    } # End Method
  } # End Buffer Loop

} # End Species Loop

# Save
final_res <- bind_rows(results_list)
write.csv(final_res, file.path(output_dir, "point_based_dual_mode_results.csv"), row.names=FALSE)
cat("Done.\n")