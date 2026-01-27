# -----------------------------------------------------------------
# Species Point-Based Experiments
# Comparison: Buffered vs Unbuffered
# -----------------------------------------------------------------

###
# 1. SETUP
###

install_now <- FALSE
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
stable_reps <- n_fit_repeats

PARAM_LOWER <- -10
PARAM_UPPER <- 10
INIT_LOWER <- -2
INIT_UPPER <- 2

res_m <- 100 
buffers_m <- c(100, 200, 500)
test_buffer_m <- 200
hex_m <- 100

# Output Directory
output_dir <- file.path("species_experiments", "output", "points")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


###
# 2. LOCAL HELPER FOR OCCU
###
# Mirrors fit_occuN_model logic but for standard occu
# Selects best model based on Negative Log Likelihood (NLL)
# Debug version of fit_occu_model
fit_occu_model <- function(umf, state_formula, obs_formula, n_reps, stable_reps, 
                           init_lower, init_upper, method = "Nelder-Mead") {
  
  n_state <- length(all.vars(state_formula)) + 1
  n_obs <- length(all.vars(obs_formula)) + 1
  n_params <- n_state + n_obs
  
  best_fm <- NULL
  min_nll <- Inf
  stable_count <- 0
  tolerance <- 0.01
  
  full_formula <- as.formula(paste(paste(deparse(obs_formula), collapse=""), 
                                   paste(deparse(state_formula), collapse="")))
  
  last_error <- "No attempts made" 
  
  for (i in 1:n_reps) {
    starts <- runif(n_params, init_lower, init_upper)
    
    # 1. Assign result to fm_try
    fm_try <- try(unmarked::occu(full_formula, data = umf, starts = starts, 
                                 method = method, se = FALSE), silent = TRUE)
    
    # 2. Check fm_try (NOT fm_rep)
    if (!inherits(fm_try, "try-error")) {
      
      curr_nll <- fm_try@negLogLike
      
      if (!is.nan(curr_nll)) {
        if (curr_nll < min_nll) {
          if (abs(min_nll - curr_nll) < tolerance) stable_count <- stable_count + 1
          else stable_count <- 0
          min_nll <- curr_nll
          best_fm <- fm_try
        } else if (abs(curr_nll - min_nll) < tolerance) {
          stable_count <- stable_count + 1
        }
      }
    } else {
      last_error <- as.character(fm_try)
    }
    if (stable_count >= stable_reps) break
  }
  
  if (is.null(best_fm)) {
    cat(paste("    [DEBUG] occu failed. Last error:", last_error, "\n"))
  }
  
  return(best_fm)
}

###
# 3. PREPROCESS RASTER
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
# 4. DATA PREPARATION (Templates)
###
cat("--- Preparing Master Data Structures (Templates) ---\n")

template_species <- species_names[1]

# A. Master Train Data (Unbuffered/Points)
train_data_res <- prepare_train_data(
  state_covs = state_cov_names,
  obs_covs = obs_cov_names,
  cov_tif = cov_tif_albers_raw,
  state_standardization_params = state_cov_params,
  placeholder_spec_name = template_species
)
master_train_df <- train_data_res$train_df


# B. Master Test Data (Unbuffered/Points)
master_test_df_unbuffered <- prepare_test_data(
  state_covs = state_cov_names,
  obs_covs = obs_cov_names,
  cov_tif = cov_tif_albers_raw, 
  standardization_params = state_cov_params,
  placeholder_spec_name = template_species
)
# C. Master Test Data (Buffered 200m - for occuN)
test_structs <- prepare_test_spatial_structures(
  test_df = master_test_df_unbuffered,
  albers_crs = albers_crs_str,
  buffer_m = test_buffer_m, 
  cov_raster_albers = cov_tif_albers,
  area_raster = area_j_raster
)
master_test_df_buffered <- test_structs$test_df
w_matrix_test <- test_structs$w_matrix

# D. Master Test Data (Buffered Mean Covs - for occu)
# D. Master Test Data (Buffered Mean Covs - for occu)
cat("--- Extracting Mean Test Covariates for Buffered occu ---\n")

# 1. Create geometries (using buffered df to ensure 'site' column exists)
test_geoms_vect <- terra::vect(voronoi_clipped_buffers(
  sf::st_as_sf(master_test_df_buffered, coords=c("longitude","latitude"), crs=4326) %>% 
    sf::st_transform(albers_crs_str), 
  buffer_dist = test_buffer_m
))

# 2. Extract Means
test_means <- terra::extract(cov_tif_albers, test_geoms_vect, fun=mean, na.rm=TRUE, ID=FALSE)

# 3. Bind Means (CRITICAL FIX: Remove original point covs first!)
master_test_df_buffered_occu <- master_test_df_buffered %>%
  dplyr::select(-all_of(state_cov_names)) %>%  # <--- DROPS DUPLICATES
  cbind(test_means)                            # <--- ADDS MEANS

###
# 5. MAIN EXPERIMENT LOOP
###

all_pred_results <- list()

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
  spec_test_obs <- get_species_obs_df(sp, 2018)
  
  # --- Update Dataframes ---
  
  # 1. Train DF (Unbuffered)
  current_train_df <- master_train_df
  current_train_df$species_observed <- NULL
  current_train_df <- inner_join(current_train_df, spec_train_obs, by="checklist_id")

  # 2. Test DF (Unbuffered)
  curr_test_unbuf <- master_test_df_unbuffered
  curr_test_unbuf$species_observed <- NULL
  curr_test_unbuf <- inner_join(curr_test_unbuf, spec_test_obs, by="checklist_id")
  
  # 3. Test DF (Buffered occuN)
  curr_test_buf <- master_test_df_buffered
  curr_test_buf$species_observed <- NULL
  curr_test_buf <- inner_join(curr_test_buf, spec_test_obs, by="checklist_id")
  
  # 4. Test DF (Buffered occu)
  curr_test_buf_occu <- master_test_df_buffered_occu
  curr_test_buf_occu$species_observed <- NULL
  curr_test_buf_occu <- inner_join(curr_test_buf_occu, spec_test_obs, by="checklist_id")
  
  
  # --- SYNCHRONIZE TEST SPLITS ---
  cat("  Generating Synchronized Test Splits...\n")
  
  test_splits_unbuf <- list()
  test_splits_buf <- list()
  test_splits_buf_occu <- list()
  
  for (r in 1:n_test_repeats) {
    # Sample IDs
    sampled_unbuf <- spatial_subsample_dataset(curr_test_unbuf, hex_m/1000, r)
    selected_ids <- sampled_unbuf$checklist_id
    
    # Store Splits
    test_splits_unbuf[[r]] <- sampled_unbuf
    test_splits_buf[[r]] <- curr_test_buf[curr_test_buf$checklist_id %in% selected_ids, ]
    test_splits_buf_occu[[r]] <- curr_test_buf_occu[curr_test_buf_occu$checklist_id %in% selected_ids, ]
  }
  
  # Define Formulas
  state_form <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))
  obs_form <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
  
  
  # =========================================================
  # EXPERIMENT A: UNBUFFERED (AAAI Style)
  # Buffer = 0 (Points only) | Model = occu ONLY
  # =========================================================
  cat("  [Experiment A] Unbuffered (Point-based, Buffer=0)\n")
  
  for (method in method_names) {
    cat(sprintf("    - Method: %s... ", method))
    
    # Filter Data
    filtered_train <- current_train_df
    if (method == "lat-long") {
      filtered_train$site <- filtered_train$locality_id
    } else if (method == "1to10") {
      filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() <= 10) %>% ungroup()
      filtered_train$site <- filtered_train$locality_id
    } else if (method == "2to10") {
      filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() >= 2, n() <= 10) %>% ungroup()
      filtered_train$site <- filtered_train$locality_id
    }

    # Prepare Data
    occu_df <- filtered_train %>%
      dplyr::select(site, species_observed, all_of(state_cov_names), all_of(obs_cov_names)) %>%
      group_by(site) %>% mutate(visit = row_number()) %>% ungroup()
    
    y_wide <- occu_df %>% dplyr::select(site, visit, species_observed) %>% 
      pivot_wider(names_from=visit, values_from=species_observed) %>% dplyr::select(-site) %>% as.matrix()
    

    site_covs <- occu_df %>% group_by(site) %>% slice(1) %>% 
      ungroup() %>%
      dplyr::select(all_of(state_cov_names)) %>% as.data.frame()
    
    obs_covs_list <- list()
    for(v in obs_cov_names){
      obs_covs_list[[v]] <- occu_df %>% dplyr::select(site, visit, all_of(v)) %>% 
        pivot_wider(names_from=visit, values_from=all_of(v)) %>% dplyr::select(-site) %>% as.matrix()
    }
    
    umf <- unmarkedFrameOccu(y = y_wide, siteCovs = site_covs, obsCovs = obs_covs_list)
    
    # FIT OCCU (With Best-Fit Logic)
    fm <- fit_occu_model(umf, state_form, obs_form, 
                         n_reps = n_fit_repeats, stable_reps = stable_reps, 
                         init_lower = INIT_LOWER, init_upper = INIT_UPPER)
    
    if (!is.null(fm)) {
      auc_vec <- numeric(n_test_repeats)
      auprc_vec <- numeric(n_test_repeats)
      
      for(r in 1:n_test_repeats){
        t_df <- test_splits_unbuf[[r]]
        X_state <- model.matrix(state_form, t_df)
        X_obs <- model.matrix(obs_form, t_df)
        
        pred_psi <- plogis(X_state %*% coef(fm, "state"))
        pred_det <- plogis(X_obs %*% coef(fm, "det"))
        
        # USE METRICS FUNCTION
        met <- calculate_classification_metrics(pred_psi * pred_det, t_df$species_observed)
        auc_vec[r] <- met$auc
        auprc_vec[r] <- met$auprc
      }
      
      # SAVE ALL REPEATS
      all_pred_results[[length(all_pred_results)+1]] <- data.frame(
        species = sp, buffer = 0, method = method, model = "occu",
        test_repeat = 1:n_test_repeats,
        auc = auc_vec,
        auprc = auprc_vec
      )
    } else {
        cat("FAILED.\n")
    }
    cat("Done.\n")
  } # End Method Loop (Unbuffered)
  
  
  # =========================================================
  # EXPERIMENT B: BUFFERED (occuN Style)
  # Buffer = 100, 200, 500
  # =========================================================
  
  for (buf in buffers_m) {
    cat(sprintf("  [Experiment B] Buffered (Size: %dm)\n", buf))
    
    for (method in method_names) {
      cat(sprintf("    - Method: %s... ", method))
      
      # Filter
      filtered_train <- current_train_df
      if (method == "lat-long") {
        filtered_train$site <- filtered_train$locality_id
      } else if (method == "1to10") {
        filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() <= 10) %>% ungroup()
        filtered_train$site <- filtered_train$locality_id
      } else if (method == "2to10") {
        filtered_train <- filtered_train %>% group_by(locality_id) %>% filter(n() >= 2, n() <= 10) %>% ungroup()
        filtered_train$site <- filtered_train$locality_id
      }
      
      # Create Geometries & W Matrix
      train_geoms <- create_site_geometries(filtered_train, cov_tif_albers, buffer_m = buf)
      train_w <- generate_overlap_matrix(train_geoms, cov_tif_albers)
      
      # Extract MEAN covariates for occu
      train_means_vect <- terra::extract(cov_tif_albers, terra::vect(train_geoms), fun=mean, na.rm=TRUE, ID=FALSE)
      train_means <- cbind(train_geoms[,c("site"), drop=F], train_means_vect)
      
      # Join Means back to filtered_train for occu
      train_df_occu_buf <- filtered_train %>%
        dplyr::select(-all_of(state_cov_names)) %>% 
        inner_join(sf::st_drop_geometry(train_means), by="site")
      
      
      # --- MODEL 1: occuN ---
      site_lookup <- filtered_train %>% dplyr::select(checklist_id, site)
      umf_n <- prepare_occuN_data(filtered_train, site_lookup, train_w, 
                                obs_cov_names, full_raster_covs)
      
      fm_n <- fit_occuN_model(umf_n, state_form, obs_form, 
                            n_reps=n_fit_repeats, stable_reps=stable_reps,
                            lower=PARAM_LOWER, upper=PARAM_UPPER,
                            init_lower=INIT_LOWER, init_upper=INIT_UPPER,
                            optimizer=selected_optimizer)
      
      if(!is.null(fm_n)) {
        auc_vec <- numeric(n_test_repeats)
        auprc_vec <- numeric(n_test_repeats)
        alphas <- coef(fm_n, 'det'); betas <- coef(fm_n, 'state')
        
        for(r in 1:n_test_repeats){
          t_df <- test_splits_buf[[r]]
          X_state <- model.matrix(state_form, t_df)
          X_obs <- model.matrix(obs_form, t_df)
          
          lambda <- exp(X_state %*% betas)
          psi <- 1 - exp(-lambda * t_df$area_j)
          p <- plogis(X_obs %*% alphas)
          
          # USE METRICS FUNCTION
          met <- calculate_classification_metrics(psi * p, t_df$species_observed)
          auc_vec[r] <- met$auc
          auprc_vec[r] <- met$auprc
        }
        
        all_pred_results[[length(all_pred_results)+1]] <- data.frame(
          species = sp, buffer = buf, method = method, model = "occuN",
          test_repeat = 1:n_test_repeats,
          auc = auc_vec,
          auprc = auprc_vec
        )
      } else {
         cat("occuN FAILED. ")
      }
      
      
      # --- MODEL 2: occu (Buffered) ---
      occu_df <- train_df_occu_buf %>%
        group_by(site) %>% mutate(visit = row_number()) %>% ungroup()
      
      y_wide <- occu_df %>% dplyr::select(site, visit, species_observed) %>% 
        pivot_wider(names_from=visit, values_from=species_observed) %>% dplyr::select(-site) %>% as.matrix()
      
      site_covs <- occu_df %>% group_by(site) %>% slice(1) %>% 
        ungroup() %>%
        dplyr::select(all_of(state_cov_names)) %>% as.data.frame()
      
      obs_covs_list <- list()
      for(v in obs_cov_names){
        obs_covs_list[[v]] <- occu_df %>% dplyr::select(site, visit, all_of(v)) %>% 
          pivot_wider(names_from=visit, values_from=all_of(v)) %>% dplyr::select(-site) %>% as.matrix()
      }
      
      umf_o <- unmarkedFrameOccu(y = y_wide, siteCovs = site_covs, obsCovs = obs_covs_list)
      
      # FIT OCCU (With Best-Fit Logic)
      fm_o <- fit_occu_model(umf_o, state_form, obs_form, 
                             n_reps = n_fit_repeats, stable_reps = stable_reps, 
                             init_lower = INIT_LOWER, init_upper = INIT_UPPER)
      
      if (!is.null(fm_o)) {
        auc_vec <- numeric(n_test_repeats)
        auprc_vec <- numeric(n_test_repeats)
        
        for(r in 1:n_test_repeats){
          t_df <- test_splits_buf_occu[[r]]
          X_state <- model.matrix(state_form, t_df)
          X_obs <- model.matrix(obs_form, t_df)
          
          pred_psi <- plogis(X_state %*% coef(fm_o, "state"))
          pred_det <- plogis(X_obs %*% coef(fm_o, "det"))
          
          # USE METRICS FUNCTION
          met <- calculate_classification_metrics(pred_psi * pred_det, t_df$species_observed)
          auc_vec[r] <- met$auc
          auprc_vec[r] <- met$auprc
        }
        
        all_pred_results[[length(all_pred_results)+1]] <- data.frame(
          species = sp, buffer = buf, method = method, model = "occu",
          test_repeat = 1:n_test_repeats,
          auc = auc_vec,
          auprc = auprc_vec
        )
      } else {
         cat("occu FAILED.")
      }
      
      cat("Done.\n")
    } # End Method
  } # End Buffer Loop
  
} # End Species Loop

# Save
final_res <- bind_rows(all_pred_results)
write.csv(final_res, file.path(output_dir, "points_results.csv"), row.names=FALSE)
cat("Done.\n")