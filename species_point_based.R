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

# Output Directory
output_dir <- file.path("species_experiments", "point_based_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

###
# 2. PREPROCESS RASTER
###

cat("--- Pre-processing raster data... ---\n")

# 1. Define CRS
albers_crs_str <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# 2. Load Native Raster
state_cov_raster_raw <- terra::rast(file.path("state_covariate_raster", "state_covariates.tif"))
terra::crs(state_cov_raster_raw) <- "+proj=longlat +datum=WGS84"
names(state_cov_raster_raw) <- state_cov_names

# 3. Project to Albers (Raw Values)
cov_tif_albers_raw <- terra::project(state_cov_raster_raw, albers_crs_str, method="bilinear", res = res_m)

# 4. STANDARDIZE RASTER
standardization_results <- standardize_state_covs(cov_tif_albers_raw)
cov_tif_albers <- standardization_results$raster
state_cov_params <- standardization_results$params

names(cov_tif_albers) <- names(cov_tif_albers_raw)

# 5. Area Raster
cell_area_km2 <- (res_m / 1000) * (res_m / 1000)
area_j_raster <- cov_tif_albers[[1]] * 0 + cell_area_km2
names(area_j_raster) <- "area"

# 6. Full Raster Covs (for occuN preparation)
full_raster_covs <- as.data.frame(terra::values(cov_tif_albers))[, state_cov_names, drop = FALSE]
full_raster_covs[is.na(full_raster_covs)] <- 0

###
# 3. PREPARE TEST STRUCTURES (Dual Mode)
###
cat("--- Preparing Master Test Data Structures ---\n")

template_species <- species_names[1]

# A. Unbuffered Test Data (Points)
# Used for Experiment A (AAAI style)
master_test_df_unbuffered <- prepare_test_data(
  state_covs = state_cov_names,
  obs_covs = obs_cov_names,
  cov_tif = cov_tif_albers_raw, 
  standardization_params = state_cov_params,
  placeholder_spec_name = template_species
)

# B. Buffered Test Data (Fixed 200m)
# Used for Experiment B (occuN style)
test_structs <- prepare_test_spatial_structures(
  test_df = master_test_df_unbuffered,
  albers_crs = albers_crs_str,
  buffer_m = test_buffer_m, 
  cov_raster_albers = cov_tif_albers,
  area_raster = area_j_raster
)
master_test_df_buffered <- test_structs$test_df
w_matrix_test <- test_structs$w_matrix

# C. Buffered Mean Covariates (for occu on buffered data)
# We need to extract the MEAN of the 200m buffer for the occu model comparison
cat("--- Extracting Mean Test Covariates for Buffered occu ---\n")
test_geoms_vect <- terra::vect(voronoi_clipped_buffers(
  sf::st_as_sf(master_test_df_unbuffered, coords=c("longitude","latitude"), crs=4326) %>% sf::st_transform(albers_crs_str), 
  buffer_dist = test_buffer_m
))
# Extract mean values, ensure column names don't clash or are handled
test_means <- terra::extract(cov_tif_albers, test_geoms_vect, fun=mean, na.rm=TRUE, ID=FALSE)
master_test_df_buffered_occu <- cbind(master_test_df_buffered, test_means) 


###
# 4. MAIN EXPERIMENT LOOP
###

results_list <- list()

# Helper to load simple observation data
get_species_obs_df <- function(sp_name, year) {
  filename <- paste0(sp_name, "_zf_filtered_region_", year, ".csv")
  fpath <- file.path("checklist_data", "species", sp_name, filename)
  df <- read.delim(fpath, sep = ",", header = TRUE)
  df <- df[!is.na(df$duration_minutes),]
  if (year == 2017) {
    df <- df[df$observation_date >= "2017-05-15" & df$observation_date <= "2017-07-09",]
  } else {
    df <- df[df$observation_date >= "2018-05-15" & df$observation_date <= "2018-07-09",]
  }
  return(df[, c("checklist_id", "species_observed")])
}

cat("\n###############################################\n")
cat("STARTING POINT-BASED MODEL FITTING LOOP\n")
cat("###############################################\n")

for (sp in species_names) {
  cat(sprintf("\n>>> PROCESSING SPECIES: %s\n", sp))
  
  # --- Load Observations ---
  spec_train_obs <- get_species_obs_df(sp, 2017)
  spec_test_obs  <- get_species_obs_df(sp, 2018)
  
  # --- Update Test DFs with current species obs ---
  
  # 1. Unbuffered Test
  curr_test_unbuf <- master_test_df_unbuffered
  curr_test_unbuf$species_observed <- NULL
  curr_test_unbuf <- inner_join(curr_test_unbuf, spec_test_obs, by="checklist_id")
  
  # 2. Buffered Test (occuN)
  curr_test_buf <- master_test_df_buffered
  curr_test_buf$species_observed <- NULL
  curr_test_buf <- inner_join(curr_test_buf, spec_test_obs, by="checklist_id")
  
  # 3. Buffered Test (occu - Mean Covs)
  curr_test_buf_occu <- master_test_df_buffered_occu
  curr_test_buf_occu$species_observed <- NULL
  curr_test_buf_occu <- inner_join(curr_test_buf_occu, spec_test_obs, by="checklist_id")
  
  
  # =========================================================
  # EXPERIMENT A: UNBUFFERED (AAAI Style)
  # Buffer = 0 (Points only) | Model = occu ONLY
  # =========================================================
  cat("  [Experiment A] Unbuffered (Point-based, Buffer=0)\n")
  
  # Prepare Base Train (Point Extraction)
  # We use prepare_train_data to get the DF structure, then overwrite observations
  train_data_res <- prepare_train_data(
    state_covs = state_cov_names,
    obs_covs = obs_cov_names,
    cov_tif = cov_tif_albers_raw,
    state_standardization_params = state_cov_params,
    placeholder_spec_name = sp # Loads specific file internally
  )
  base_train_unbuf <- train_data_res$train_df
  
  # Generate Test Repeats (Subsampling) for Unbuffered
  test_splits_unbuf <- list()
  for (r in 1:n_test_repeats) {
    test_splits_unbuf[[r]] <- spatial_subsample_dataset(curr_test_unbuf, hex_m/1000, r)
  }
  
  for (method in method_names) {
    cat(sprintf("    - Method: %s... ", method))
    
    # Filter Data
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
    # Pivot for unmarked (Site x Visit)
    occu_df <- filtered_train %>%
      dplyr::select(site, species_observed, all_of(state_cov_names), all_of(obs_cov_names)) %>%
      group_by(site) %>% mutate(visit = row_number()) %>% ungroup()
    
    y_wide <- occu_df %>% dplyr::select(site, visit, species_observed) %>% 
      pivot_wider(names_from=visit, values_from=species_observed) %>% dplyr::select(-site) %>% as.matrix()
    
    site_covs <- occu_df %>% group_by(site) %>% slice(1) %>% 
      dplyr::select(all_of(state_cov_names)) %>% as.data.frame()
    
    obs_covs_list <- list()
    for(v in obs_cov_names){
      obs_covs_list[[v]] <- occu_df %>% dplyr::select(site, visit, all_of(v)) %>% 
        pivot_wider(names_from=visit, values_from=all_of(v)) %>% dplyr::select(-site) %>% as.matrix()
    }
    
    umf <- unmarkedFrameOccu(y = y_wide, siteCovs = site_covs, obsCovs = obs_covs_list)
    
    # Fit
    state_form <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))
    obs_form   <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
    
    fm <- try(unmarked::occu(as.formula(paste(paste(deparse(obs_form), collapse=""), paste(deparse(state_form), collapse=""))), umf), silent = TRUE)
    
    if (!inherits(fm, "try-error")) {
      # Predict on Unbuffered Test Data
      auc_vec <- c(); auprc_vec <- c()
      
      for(r in 1:n_test_repeats){
        t_df <- test_splits_unbuf[[r]]
        X_state <- model.matrix(state_form, t_df)
        X_obs   <- model.matrix(obs_form, t_df)
        
        pred_psi <- plogis(X_state %*% coef(fm, "state"))
        pred_det <- plogis(X_obs %*% coef(fm, "det"))
        
        met <- calculate_classification_metrics(pred_psi * pred_det, t_df$species_observed)
        auc_vec <- c(auc_vec, met$auc)
        auprc_vec <- c(auprc_vec, met$auprc)
      }
      
      results_list[[length(results_list)+1]] <- data.frame(
        Species = sp, Buffer = 0, Method = method, Model = "occu",
        AUC = mean(auc_vec, na.rm=TRUE), AUPRC = mean(auprc_vec, na.rm=TRUE)
      )
    }
    cat("Done.\n")
  } # End Method Loop (Unbuffered)
  
  
  # =========================================================
  # EXPERIMENT B: BUFFERED (occuN Style)
  # Buffer = 100, 200, 500
  # =========================================================
  
  # Generate Test Repeats for Buffered Sets
  test_splits_buf <- list()
  test_splits_buf_occu <- list()
  
  for (r in 1:n_test_repeats) {
    # We rely on spatial_subsample_dataset using checklist_id or similar to keep rows consistent
    # But since these DFs differ in columns, we subsample them independently but using same seed logic if needed.
    # Actually spatial_subsample_dataset uses hex binning.
    test_splits_buf[[r]]      <- spatial_subsample_dataset(curr_test_buf, hex_m/1000, r)
    test_splits_buf_occu[[r]] <- spatial_subsample_dataset(curr_test_buf_occu, hex_m/1000, r)
  }
  
  for (buf in buffers_m) {
    cat(sprintf("  [Experiment B] Buffered (Size: %dm)\n", buf))
    
    # Reload Base Train to ensure clean start
    base_train_raw <- train_data_res$train_df
    
    for (method in method_names) {
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
      # Using standard helper logic but dynamically with 'buf'
      train_geoms <- create_site_geometries(filtered_train, cov_tif_albers, buffer_m = buf)
      train_w <- generate_overlap_matrix(train_geoms, cov_tif_albers)
      
      # Extract MEAN covariates for occu
      train_means_vect <- terra::extract(cov_tif_albers, terra::vect(train_geoms), fun=mean, na.rm=TRUE, ID=FALSE)
      train_means <- cbind(train_geoms[,c("site"), drop=F], train_means_vect)
      
      # Join Means back to filtered_train for occu
      train_df_occu_buf <- filtered_train %>%
        dplyr::select(-all_of(state_cov_names)) %>% # Drop point covs
        inner_join(sf::st_drop_geometry(train_means), by="site") # Add mean covs
      
      
      # --- MODEL 1: occuN ---
      tryCatch({
        site_lookup <- filtered_train %>% dplyr::select(checklist_id, site)
        umf <- prepare_occuN_data(filtered_train, site_lookup, train_w, 
                                  obs_cov_names, full_raster_covs)
        
        state_form <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))
        obs_form   <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
        
        fm <- fit_occuN_model(umf, state_form, obs_form, 
                              n_reps=n_fit_repeats, stable_reps=n_fit_repeats,
                              optimizer=selected_optimizer)
        
        if(!is.null(fm)) {
          auc_vec <- c(); auprc_vec <- c()
          
          alphas <- coef(fm, 'det'); betas  <- coef(fm, 'state')
          
          for(r in 1:n_test_repeats){
            t_df <- test_splits_buf[[r]]
            X_state <- model.matrix(state_form, t_df)
            X_obs   <- model.matrix(obs_form, t_df)
            
            lambda <- exp(X_state %*% betas)
            psi    <- 1 - exp(-lambda * t_df$area_j) # AREA DEPENDENT
            p      <- plogis(X_obs %*% alphas)
            
            met <- calculate_classification_metrics(psi * p, t_df$species_observed)
            auc_vec <- c(auc_vec, met$auc); auprc_vec <- c(auprc_vec, met$auprc)
          }
          
          results_list[[length(results_list)+1]] <- data.frame(
            Species = sp, Buffer = buf, Method = method, Model = "occuN",
            AUC = mean(auc_vec, na.rm=TRUE), AUPRC = mean(auprc_vec, na.rm=TRUE)
          )
        }
      }, error = function(e) cat("occuN err\n"))
      
      
      # --- MODEL 2: occu (Buffered) ---
      tryCatch({
        # Use train_df_occu_buf (Mean Covs)
        occu_df <- train_df_occu_buf %>%
          group_by(site) %>% mutate(visit = row_number()) %>% ungroup()
        
        y_wide <- occu_df %>% dplyr::select(site, visit, species_observed) %>% 
          pivot_wider(names_from=visit, values_from=species_observed) %>% dplyr::select(-site) %>% as.matrix()
        
        site_covs <- occu_df %>% group_by(site) %>% slice(1) %>% 
          dplyr::select(all_of(state_cov_names)) %>% as.data.frame()
        
        obs_covs_list <- list()
        for(v in obs_cov_names){
          obs_covs_list[[v]] <- occu_df %>% dplyr::select(site, visit, all_of(v)) %>% 
            pivot_wider(names_from=visit, values_from=all_of(v)) %>% dplyr::select(-site) %>% as.matrix()
        }
        
        umf <- unmarkedFrameOccu(y = y_wide, siteCovs = site_covs, obsCovs = obs_covs_list)
        fm <- try(unmarked::occu(as.formula(paste(paste(deparse(obs_form), collapse=""), paste(deparse(state_form), collapse=""))), umf), silent=TRUE)
        
        if (!inherits(fm, "try-error")) {
          auc_vec <- c(); auprc_vec <- c()
          
          for(r in 1:n_test_repeats){
            t_df <- test_splits_buf_occu[[r]]
            X_state <- model.matrix(state_form, t_df)
            X_obs   <- model.matrix(obs_form, t_df)
            
            pred_psi <- plogis(X_state %*% coef(fm, "state"))
            pred_det <- plogis(X_obs %*% coef(fm, "det"))
            
            met <- calculate_classification_metrics(pred_psi * pred_det, t_df$species_observed)
            auc_vec <- c(auc_vec, met$auc); auprc_vec <- c(auprc_vec, met$auprc)
          }
          
          results_list[[length(results_list)+1]] <- data.frame(
            Species = sp, Buffer = buf, Method = method, Model = "occu",
            AUC = mean(auc_vec, na.rm=TRUE), AUPRC = mean(auprc_vec, na.rm=TRUE)
          )
        }
      }, error = function(e) cat("occu err\n"))
      
      cat("Done.\n")
    } # End Method
  } # End Buffer Loop
  
} # End Species Loop

# Save
final_res <- bind_rows(results_list)
write.csv(final_res, file.path(output_dir, "point_based_results.csv"), row.names=FALSE)
cat("Done.\n")