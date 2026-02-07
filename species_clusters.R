# -----------------------------------------------------------------
# Real Species Experiments for occuN model
#
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
library(PRROC)
library(terra)

source(file.path("R", "utils.R"))
source(file.path("R", "data_helpers.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))
source(file.path("R", "plotting_helpers.R"))
# Added SLIC source
source(file.path("R", "clustering", "slic.R"))

set.seed(123) 

###
# 2. CONFIGS
###

# Species list
species_names <- c(
  "AMCR", "AMRO", "BAEA", "BKHGRO", "BRCR", "BUTI", "CASC", "CHBCHI", 
  "COHA", "HAFL", "HAWO", "HEWA", "MAWA", "MOQU", "NOFL", "NOOW", 
  "OLFL", "PAFL", "PAWR", "PIWO", "REHA", "SOSP", "SPTO", "SWTH", 
  "WAVI", "WEPE", "WETA", "WIWA", "WRENTI", "YEBCHA", "YEWA"
)

# Comparison methods
method_names <- c(
  "1to10", 
  "2to10", 
  "2to10-sameObs", 
  "1-kmSq",
  "lat-long", 
  "rounded-4", 
  "clustGeo-50-5",
  "clustGeo-50-10",
  "clustGeo-50-20",
  "clustGeo-50-30",
  "clustGeo-50-40",
  "clustGeo-50-50",
  "clustGeo-50-55",
  "clustGeo-50-60",
  "clustGeo-50-70",
  "clustGeo-50-80",
  "clustGeo-50-90",
  "clustGeo-50-95",
  "DBSC",
  "BayesOptClustGeo"
  # "SLIC-0.05-3200",
  # "SLIC-0.1-3200",
  # "SLIC-0.15-3200",
  # "SLIC-0.2-3200",
  # "SLIC-0.25-3200",
  # "SLIC-0.3-3200",
  # "SLIC-0.35-3200",
  # "SLIC-0.4-3200",
  # "SLIC-0.45-3200"
)

# Methods to plot
methods_to_plot <- c(
  "1to10", "2to10", "2to10-sameObs", "lat-long", "1-kmSq", "rounded-4", "DBSC"
)

methods_to_plot_clustGeo <- c(
  "clustGeo-50-5",
  "clustGeo-50-10",
  "clustGeo-50-20",
  "clustGeo-50-30",
  "clustGeo-50-40",
  "clustGeo-50-50",
  "clustGeo-50-55",
  "clustGeo-50-60",
  "clustGeo-50-70",
  "clustGeo-50-80",
  "clustGeo-50-90",
  "clustGeo-50-95"
)

methods_to_plot_slic <- c(
  "SLIC-0.05-3200",
  "SLIC-0.1-3200",
  "SLIC-0.15-3200",
  "SLIC-0.2-3200",
  "SLIC-0.25-3200",
  "SLIC-0.3-3200",
  "SLIC-0.35-3200",
  "SLIC-0.4-3200",
  "SLIC-0.45-3200"
)

# Covariates
state_cov_names <- c("elevation", "TCB", "TCG", "TCW", "TCA")
obs_cov_names <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")

# Optimization & Simulation Settings
selected_optimizer <- "nlminb"
n_fit_repeats <- 50
n_test_repeats <- 25

res_m <- 100 
buffer_m <- 200
hex_m <- 100

PARAM_LOWER <- -10
PARAM_UPPER <- 10
INIT_LOWER <- -2
INIT_UPPER <- 2


min_uniloc_points <- 2 # Options: 1, 2, 3, ...
max_uniloc_points <- 10 # Options: 1, 3, ... 10, ... "all"

# Output Directory
output_dir <- file.path("species_experiments", "output", "clusters")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


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

# Boundary for plotting
boundary_shapefile_path <- file.path("state_covariate_raster", "boundary", "boundary.shp")


###
# 4. VERIFY CHECKLIST CONSISTENCY
###
cat("\n###############################################\n")
cat("VERIFYING DATA CONSISTENCY ACROSS SPECIES\n")
cat("###############################################\n")

# Helper to pull checklist IDs for a species/year without doing full processing
get_checklist_ids <- function(sp_name, year) {
    filename <- paste0(sp_name, "_zf_filtered_region_", year, ".csv")
    fpath <- file.path("checklist_data", "species", sp_name, filename)
    if (!file.exists(fpath)) stop(paste("File not found:", fpath))
    
    df <- read.delim(fpath, sep = ",", header = TRUE)
    # Apply exactly the same filters as prepare_train/test_data
    df <- df[!is.na(df$duration_minutes),]
    if (year == 2017) {
        df <- df[df$observation_date >= "2017-05-15" & df$observation_date <= "2017-07-09",]
    } else {
        df <- df[df$observation_date >= "2018-05-15" & df$observation_date <= "2018-07-09",]
    }
    return(sort(df$checklist_id))
}

cat("Reading reference IDs from first species:", species_names[1], "\n")
ref_train_ids <- get_checklist_ids(species_names[1], 2017)
ref_test_ids <- get_checklist_ids(species_names[1], 2018)

cat("Checking all other species...\n")
for (sp in species_names[-1]) {
    cur_train <- get_checklist_ids(sp, 2017)
    cur_test <- get_checklist_ids(sp, 2018)
    
    if (!identical(ref_train_ids, cur_train)) stop(paste("Train Checklist IDs do not match for", sp))
    if (!identical(ref_test_ids, cur_test)) stop(paste("Test Checklist IDs do not match for", sp))
}
cat("SUCCESS: All species share identical checklist_ids for Train (2017) and Test (2018) sets.\n")


###
# 5. SPATIAL PROCESSING
###
cat("\n###############################################\n")
cat("RUNNING SPATIAL PROCESSING\n")
cat("Using template species:", species_names[1], "\n")
cat("###############################################\n")

template_species <- species_names[1]

# === 5.1 DATA PREPARATION (Template) ===
cat("--- Preparing Master Train/Test Data Structures ---\n")

# This creates the "skeleton" dataframe with correct environmental covariates and rows.
# The 'species_observed' column here will be overwritten in the loop, but used for initial clustering.
train_data_res <- prepare_train_data(
    state_covs = state_cov_names, 
    obs_covs = obs_cov_names, 
    cov_tif = cov_tif_albers_raw, 
    state_standardization_params = state_cov_params,
    placeholder_spec_name = template_species
)

master_train_df <- train_data_res$train_df
full_standardization_params <- train_data_res$standardization_params

master_test_df <- prepare_test_data(
    state_covs = state_cov_names, 
    obs_covs = obs_cov_names, 
    cov_tif = cov_tif_albers_raw, 
    standardization_params = full_standardization_params,
    placeholder_spec_name = template_species
)
# === 5.2 CLUSTERING (Modified for Option ii) ===
cat("--- Computing clusterings ---\n")

# 1. Define the sets of methods

# ONLY keep methods here that MUST see singletons (1 obs/site)
baseline_methods <- c("1to10", "lat-long")

# ALL methods that should compare the "2-to-10" filtered dataset
spatial_methods <- c(
  "2to10",              # Moved here to ensure identical points
  "2to10-sameObs",      # Moved here to ensure identical points
  "1-kmSq",             # Moved per your request
  "rounded-4",          # Moved per your request
  "DBSC", 
  "BayesOptClustGeo"
)

# Add the specific parameter versions if they exist in your 'method_names' list
spatial_methods <- c(spatial_methods, 
                     grep("clustGeo-", method_names, value = TRUE),
                     grep("SLIC-", method_names, value = TRUE))

# Intersect with what the user actually requested in 'method_names'
baseline_run_list <- intersect(method_names, baseline_methods)
spatial_run_list  <- intersect(method_names, spatial_methods)

# 2. Create the "Option (ii)" Dataset
# Logic: Keep only locations with >= 2 visits, then downsample to max 10.
cat("--- Creating Filtered Input for Spatial Clustering (Min 2, Max 10 per Loc) ---\n")

train_df_spatial_input <- master_train_df %>%
  group_by(locality_id) %>%
  filter(n() >= min_uniloc_points) %>%       # STEP A: Remove single-visit locations
  slice_sample(n = max_uniloc_points) %>%   # STEP B: Downsample to max 10 (like 2to10)
  ungroup()

# 3. Run Clustering in Two Batches

# Batch A: Baselines (Run on FULL data)
# (Note: '2to10' handles its own filtering inside auk, so we give it everything)
results_baseline <- get_clusterings(baseline_run_list, master_train_df, state_cov_names, NULL, cov_tif_albers)

# Batch B: Spatial Methods (Run on FILTERED data)
# Since inputs have >=2 obs, all resulting clusters will have >=2 obs.
results_spatial  <- get_clusterings(spatial_run_list, train_df_spatial_input, state_cov_names, NULL, cov_tif_albers)

# 4. Combine Results
all_clusterings <- c(results_baseline, results_spatial)

cat("--- Computing initial site geometries ---\n")
all_site_geometries <- list()

for (method_name in method_names) {
    # Retrieve the raw output from get_clusterings
    cluster_obj <- all_clusterings[[method_name]]
    
    # 1. Handle SLIC (Pre-calculated geometries)
    if (grepl("SLIC", method_name) && is.list(cluster_obj) && "site_geoms" %in% names(cluster_obj)) {
        
        cat(sprintf("Using pre-calculated polygons for %s...\n", method_name))
        all_site_geometries[[method_name]] <- cluster_obj$site_geoms
        
    # 2. Handle Standard Methods (Buffer points)
    } else {
        # Extract just the dataframe if it's wrapped in a list (like BayesOpt or the new SLIC structure)
        if (is.list(cluster_obj) && "result_df" %in% names(cluster_obj)) {
            cluster_df <- cluster_obj$result_df
        } else {
            cluster_df <- cluster_obj
        }
        
        if (!is.null(cluster_df)) {
            all_site_geometries[[method_name]] <- create_site_geometries(cluster_df, cov_tif_albers, buffer_m, method_name)
        }
    }
}

# Record Stats PRE Split
cat("--- Summarizing Stats PRE ---\n")
all_clustering_stats_pre <- list()
stats_pre <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
stats_pre$species <- "GLOBAL_TEMPLATE" # Marked as template
all_clustering_stats_pre[[1]] <- stats_pre

# --- PLOTTING PRE-SPLIT ---
cat("--- Plotting sites (PRE-SPLIT) ---\n")
plot_sites(
    base_train_df = master_train_df,
    all_clusterings = all_clusterings,
    all_site_geometries = all_site_geometries,
    elevation_raster = cov_tif_albers_raw, 
    methods_to_plot = methods_to_plot,
    boundary_shp_path = boundary_shapefile_path,
    output_path = file.path(output_dir, "sites_PRE.png"),
    cluster_labels = TRUE
)
plot_sites(
    base_train_df = master_train_df,
    all_clusterings = all_clusterings,
    all_site_geometries = all_site_geometries,
    elevation_raster = cov_tif_albers_raw, 
    methods_to_plot = methods_to_plot_clustGeo,
    boundary_shp_path = boundary_shapefile_path,
    output_path = file.path(output_dir, "sites_clustGeo_PRE.png"),
    cluster_labels = TRUE
)
# plot_sites(
#     base_train_df = master_train_df,
#     all_clusterings = all_clusterings,
#     all_site_geometries = all_site_geometries,
#     elevation_raster = cov_tif_albers_raw, 
#     methods_to_plot = methods_to_plot_slic,
#     boundary_shp_path = boundary_shapefile_path,
#     output_path = file.path(output_dir, "sites_SLIC_PRE.png"),
#     cluster_labels = TRUE
# )


# === 5.3 SPLIT DISJOINT SITES ===
cat("--- Splitting disjoint geometries ---\n")
for (method_name in names(all_site_geometries)) {
    curr_geoms <- all_site_geometries[[method_name]]
    curr_data_obj <- all_clusterings[[method_name]]
    
    is_list_obj <- is.list(curr_data_obj) && "result_df" %in% names(curr_data_obj)
    curr_data <- if(is_list_obj) curr_data_obj$result_df else curr_data_obj
    
    split_res <- disjoint_site_geometries(curr_geoms, curr_data)
    all_site_geometries[[method_name]] <- split_res$geoms
    
    if (is_list_obj) {
        all_clusterings[[method_name]]$result_df <- split_res$data
    } else {
        all_clusterings[[method_name]] <- split_res$data
    }
}

# Record Stats POST Split
cat("--- Summarizing Stats POST ---\n")
all_clustering_stats_post <- list()
stats_post <- summarize_clusterings(all_clusterings, all_site_geometries, units = "km")
stats_post$species <- "GLOBAL_TEMPLATE"
all_clustering_stats_post[[1]] <- stats_post

# --- PLOTTING POST-SPLIT (Run Once) ---
cat("--- Plotting sites (POST-SPLIT) ---\n")
plot_sites(
    base_train_df = master_train_df,
    all_clusterings = all_clusterings,
    all_site_geometries = all_site_geometries,
    elevation_raster = cov_tif_albers_raw, 
    methods_to_plot = methods_to_plot,
    boundary_shp_path = boundary_shapefile_path,
    output_path = file.path(output_dir, "sites_POST.png"),
    cluster_labels = TRUE
)
plot_sites(
    base_train_df = master_train_df,
    all_clusterings = all_clusterings,
    all_site_geometries = all_site_geometries,
    elevation_raster = cov_tif_albers_raw, 
    methods_to_plot = methods_to_plot_clustGeo,
    boundary_shp_path = boundary_shapefile_path,
    output_path = file.path(output_dir, "sites_clustGeo_POST.png"),
    cluster_labels = TRUE
)

# plot_sites(
#     base_train_df = master_train_df,
#     all_clusterings = all_clusterings,
#     all_site_geometries = all_site_geometries,
#     elevation_raster = cov_tif_albers_raw, 
#     methods_to_plot = methods_to_plot_slic,
#     boundary_shp_path = boundary_shapefile_path,
#     output_path = file.path(output_dir, "sites_SLIC_POST.png"),
#     cluster_labels = TRUE
# )

# === 5.4 W MATRICES ===
cat("--- Generating W matrices ---\n")
all_w_matrices <- list()
for (m_name in names(all_site_geometries)) {
    if (!is.null(all_site_geometries[[m_name]])) {
        all_w_matrices[[m_name]] <- generate_overlap_matrix(all_site_geometries[[m_name]], cov_tif_albers)
    }
}

# === 5.5 TEST SPATIAL STRUCTURES ===
cat("--- Preparing Test Structures ---\n")
test_structures <- prepare_test_spatial_structures(
    test_df = master_test_df,
    albers_crs = albers_crs_str,
    buffer_m = buffer_m,
    cov_raster_albers = cov_tif_albers,
    area_raster = area_j_raster
)
master_test_df_ready <- test_structures$test_df
w_matrix_test <- test_structures$w_matrix


###
# 6. MODEL FITTING LOOP (Per Species)
###

all_param_results <- list()
all_pred_results <- list()

# Helper to load simple observation data
get_species_obs_df <- function(sp_name, year) {
    filename <- paste0(sp_name, "_zf_filtered_region_", year, ".csv")
    fpath <- file.path("checklist_data", "species", sp_name, filename)
    df <- read.delim(fpath, sep = ",", header = TRUE)
    # Apply filters
    df <- df[!is.na(df$duration_minutes),]
    if (year == 2017) {
        df <- df[df$observation_date >= "2017-05-15" & df$observation_date <= "2017-07-09",]
    } else {
        df <- df[df$observation_date >= "2018-05-15" & df$observation_date <= "2018-07-09",]
    }
    # Return minimal DF for joining
    return(df[, c("checklist_id", "species_observed")])
}

cat("\n###############################################\n")
cat("STARTING SPECIES MODEL FITTING LOOP\n")
cat("###############################################\n")

for (species_name in species_names) {
  
  cat(sprintf("\n>>> PROCESSING SPECIES: %s\n", species_name))

  # --- Update Train Data ---
  # Load fresh obs for this species
  spec_train_obs <- get_species_obs_df(species_name, 2017)
  
  # Update master_train_df by joining on checklist_id
  # We use inner_join to be safe, but we know IDs match. 
  # We select only the columns we need to replace to avoid duplicate column names.
  current_train_df <- master_train_df
  
  # Remove the old 'species_observed' (from template or previous loop)
  current_train_df$species_observed <- NULL
  
  # Join new obs
  current_train_df <- inner_join(current_train_df, spec_train_obs, by = "checklist_id")
  
  
  # --- Update Test Data ---
  spec_test_obs <- get_species_obs_df(species_name, 2018)
  current_test_df <- master_test_df_ready
  current_test_df$species_observed <- NULL
  current_test_df <- inner_join(current_test_df, spec_test_obs, by = "checklist_id")
  
  
  cat(sprintf("Training Points: %d; Test Points: %d\n", nrow(spec_train_obs), nrow(spec_test_obs)))

  # --- Regenerate Test Splits (Fast) ---
  # Test splits are just subsamples of rows. We must regenerate them because
  # the 'species_observed' column in the source (current_test_df) has changed.
  test_splits_list <- list()
  for (r in 1:n_test_repeats) {
    test_splits_list[[r]] <- spatial_subsample_dataset(current_test_df, hex_m/1000, r)
  }

  
  # --- FIT MODELS ---
  for (method_name in method_names) {
    cat(sprintf("  - Method: %s... ", method_name))
    
    # 1. Retrieve Pre-Calculated Clustering Data
    # The clustering (site assignments) comes from the global run.
    # However, 'current_clustering_df' is derived from training data.
    # In 'get_clusterings', the output often contains the original data + 'site' column.
    # We need to ensure that the 'species_observed' in the clustering DF matches the current species.
    
    global_clust_obj <- all_clusterings[[method_name]]
    
    # Extract the DF
    if (is.list(global_clust_obj) && "result_df" %in% names(global_clust_obj)) {
      clust_df <- global_clust_obj$result_df
    } else {
      clust_df <- global_clust_obj
    }
    
    if (is.null(clust_df)) {
        cat("Skipping (NULL data)\n")
        next
    }
    
    # CRITICAL: The 'clust_df' has the template species 'species_observed'.
    # We must update it with the current species observations.
    # We can do this by joining 'spec_train_obs' again or transferring from current_train_df
    clust_df$species_observed <- NULL
    clust_df <- inner_join(clust_df, spec_train_obs, by = "checklist_id")



    if (method_name == "BayesOptClustGeo") {
        cat("  Optimizing specific clusters (GLM Proxy)... ")
        
        # 1. Prepare Data Splits
        # 'train_df_spatial_input' (min 2, max 10) is already defined globally,
        # BUT we need to make sure it has the CURRENT species observations attached.
        
        # Re-attach current species obs to the spatial input skeleton
        current_spatial_input <- train_df_spatial_input %>%
            select(-any_of("species_observed")) %>%
            inner_join(spec_train_obs, by = "checklist_id")

        # Create Validation Set (Singletons + Overflow)
        # Everything in current_train_df that is NOT in the spatial input
        validation_df <- current_train_df %>%
            filter(!checklist_id %in% current_spatial_input$checklist_id)

        # 2. Run Optimization
        opt_res <- bayesianOptimizedClustGeo(
            train_data = current_spatial_input,
            validation_data = validation_df,
            state_covs = state_cov_names,
            obs_covs = obs_cov_names,
            n_init = 10,       # Initial random points
            n_iter = 15,       # Optimization steps (adjust for time)
            env_weight = 5
        )
        
        best_rho <- opt_res$Best_Pars$rho
        best_kappa <- opt_res$Best_Pars$kappa
        cat(sprintf("[Rho: %.2f, Kappa: %.2f] ", best_rho, best_kappa))
        
        # 3. Generate Final Clusters (Re-run one last time with best params)
        # Get unique locations
        uniq_loc_df <- dplyr::distinct(current_spatial_input, latitude, longitude, .keep_all = TRUE)
        # normalize and dist (same as inside the function, simplified here for generation)
        # ... actually, relying on 'clustGeoSites' helper is cleaner if it handles raw data:
        
        # We need to manually generate the site IDs on the unique locations
        # to ensure consistency with how BayesOpt did it.
        # It's safest to just rely on clustGeoSites if it accepts raw inputs:
        final_sites_df <- clustGeoSites(
            alpha = best_rho,
            checklists = uniq_loc_df,
            state_covs = state_cov_names,
            num_sites = max(2, round(nrow(uniq_loc_df) * (best_kappa / 100.0)))
        )
        
        # 4. Map back to full spatial input
        site_lookup <- final_sites_df[, c("locality_id", "site")]
        
        # This 'clust_df' is what occuN expects
        # We join it to 'current_spatial_input'
        clust_df <- left_join(current_spatial_input, site_lookup, by = "locality_id")
        
        # Update 'global_clust_obj' so the standard code below uses it
        global_clust_obj <- clust_df 
        
        # 5. Generate W Matrix specifically for this new clustering
        # (We must do this here because the sites didn't exist before)
        # Extract geometries
        current_geoms <- create_site_geometries(clust_df, cov_tif_albers, buffer_m, "BayesOpt_Temp")
        # Handle disjoint
        split_res <- disjoint_site_geometries(current_geoms, clust_df)
        current_geoms <- split_res$geoms
        clust_df <- split_res$data # Update data with split sites
        
        # Generate W
        w_matrix <- generate_overlap_matrix(current_geoms, cov_tif_albers)
        
        # Override the standard loop variable
        # (The loop usually expects 'all_w_matrices[[method_name]]')
        # We just set w_matrix directly here, and the code below uses 'w_matrix'
        
    } else {
        # [STANDARD LOGIC] for all other methods
        global_clust_obj <- all_clusterings[[method_name]]
        
        # ... (rest of standard logic to extract clust_df) ...
        if (is.list(global_clust_obj) && "result_df" %in% names(global_clust_obj)) {
          clust_df <- global_clust_obj$result_df
        } else {
          clust_df <- global_clust_obj
        }
        
        # Ensure species obs are updated
        if(!is.null(clust_df)) {
             clust_df$species_observed <- NULL
             clust_df <- inner_join(clust_df, spec_train_obs, by = "checklist_id")
        }

        w_matrix <- all_w_matrices[[method_name]]
    }

    # # 2. Retrieve W Matrix
    # w_matrix <- all_w_matrices[[method_name]]
    
    if (is.null(w_matrix)) {
       cat("Skipping (No W matrix)\n"); next
    }
    
    # 3. Prepare UMF
    # Note: prepare_occuN_data uses 'current_train_df' (base) and 'clust_df' (sites).
    # Both now have the correct species_observed.
    umf <- prepare_occuN_data(current_train_df, clust_df, w_matrix, obs_cov_names, full_raster_covs)
    
    obs_formula <- as.formula(paste("~", paste(obs_cov_names, collapse = " + ")))
    state_formula <- as.formula(paste("~", paste(state_cov_names, collapse = " + ")))
    
    # 4. Fit Model
    fm <- fit_occuN_model(
      umf, state_formula, obs_formula,
      n_reps = n_fit_repeats, stable_reps = n_fit_repeats,
      optimizer = selected_optimizer, lower = PARAM_LOWER, upper = PARAM_UPPER,
      init_lower = INIT_LOWER, init_upper = INIT_UPPER
    )
    
    if (is.null(fm)) {
      cat("FAILED.\n")
      next
    }
    
    # 5. Save Parameters
    est_alphas <- coef(fm, 'det')   
    est_betas <- coef(fm, 'state')  
    
    state_col_names <- c("state_intercept", state_cov_names)
    obs_col_names   <- c("obs_intercept", obs_cov_names)
    
    param_row <- data.frame(
      species = species_name,
      method = method_name,
      nll = fm@negLogLike,
      convergence = 0
    )
    
    for(i in seq_along(state_col_names)) param_row[[state_col_names[i]]] <- est_betas[i]
    for(i in seq_along(obs_col_names)) param_row[[obs_col_names[i]]] <- est_alphas[i]
    
    all_param_results[[length(all_param_results) + 1]] <- param_row
    
    
    # 6. Predict & Test
    # (Uses the pre-generated test_splits_list which now contains current species data)
    auc_list <- c()
    auprc_list <- c()
    
    for (repeat_num in 1:n_test_repeats) {
      test_df <- test_splits_list[[repeat_num]]
      
      X_state <- model.matrix(state_formula, data = test_df)
      X_obs <- model.matrix(obs_formula, data = test_df)
      
      pred_psi <- 1 - exp(-(exp(X_state %*% est_betas) * test_df$area_j))
      pred_det <- plogis(X_obs %*% est_alphas)
      pred_obs_prob <- pred_psi * pred_det 
      
      metrics <- calculate_classification_metrics(pred_obs_prob, test_df$species_observed)
      
      all_pred_results[[length(all_pred_results) + 1]] <- data.frame(
        species = species_name,
        method = method_name,
        test_repeat = repeat_num,
        auc = metrics$auc,
        auprc = metrics$auprc
      )
    }
    
    cat("Done.\n")
    
    # Memory Cleanup (per method)
    rm(umf, fm)
    gc()
    
  } # End Method Loop
  
  # Memory Cleanup (per species)
  rm(test_splits_list, current_train_df, current_test_df, spec_train_obs, spec_test_obs)
  gc()
  
} # End Species Loop

# Save Final Results
write.csv(dplyr::bind_rows(all_param_results), file.path(output_dir, "estimated_parameters.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_pred_results), file.path(output_dir, "predictive_performance.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_clustering_stats_pre), file.path(output_dir, "clustering_stats_PRE.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_clustering_stats_post), file.path(output_dir, "clustering_stats_POST.csv"), row.names = FALSE)

cat("Done.\n")