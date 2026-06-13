################################################################
# Updated Simulation Experiments: Clustering Comparison
# Testing MAUP via "Cell-as-the-Micro-Truth" approach
################################################################

###
# 1. SETUP & LIBRARIES
###

install_now = FALSE
if (install_now){
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  if (!requireNamespace("devtools", quietly = FALSE)) install.packages("devtools")
  suppressMessages(devtools::install_github("anonymous97331/unmarked", ref = "main", force = TRUE))
}

library(unmarked)
library(ggplot2)
library(patchwork)
library(terra) 
library(Matrix) 
library(sf)
library(dplyr)

# Source Helpers (Ensure these point to your actual helper files)
source(file.path("R", "utils.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))

set.seed(123) 

##########
# 2. CONFIGURATION
##########

# --- Simulation repetitions ---
n_sims <- 2 # SET TO 100 FOR FULL RUN
n_reps <- 30 # Number of random-start repetitions for model fitting

# --- Landscape Parameters ---
full_grid_dim <- 200
full_n_cells <- full_grid_dim * full_grid_dim
res_m <- 100 # Each cell is 100m x 100m
cell_area_km2 <- (res_m / 1000)^2

# --- Checklist Simulation Parameters ---
total_checklists <- 5000
n_unique_locations <- 1500 

# --- True parameter values ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
true_betas <- c(beta_int = -5.0, beta_cov = 1.0)

# --- Methods to Test ---
methods_to_test <- c(
  "lat-long", "1to10", "2to10", 
  "0.25-kmSq", "1-kmSq", "4-kmSq",  # <--- EXACT MATCH FOR 5, 10, 20 CELLS
  "clustGeo-50-50", "DBSC" 
)

buffer_m <- 200
PARAM_LOWER <- -10
PARAM_UPPER <- 10
INIT_LOWER <- -2
INIT_UPPER <- 2

selected_optimizer <- "nlminb"

output_dir <- file.path("output", "simulation_experiments", "clustering_comparison")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


##########
# 3. HELPER: SIMULATE OBSERVER LOCATIONS
##########
simulate_ebird_checklists <- function(total_n, n_unique, max_coord) {
  
  # 1. Generate unique random coordinates in meters
  locs_x <- runif(n_unique, 5, max_coord - 5)
  locs_y <- runif(n_unique, 5, max_coord - 5)
  
  # 2. Generate skewed visitation frequency
  weights <- (1:n_unique)^(-1.5) 
  weights <- weights / sum(weights)
  sampled_loc_indices <- sample(1:n_unique, size = total_n, replace = TRUE, prob = weights)
  
  # 3. FIX: Convert meter coordinates to actual Lat/Longs (EPSG:4326) 
  # This ensures they pass safely through your model_helpers.R pipeline
  sim_coords <- data.frame(
    x = locs_x[sampled_loc_indices],
    y = locs_y[sampled_loc_indices]
  )
  sim_sf <- sf::st_as_sf(sim_coords, coords = c("x", "y"), crs = "EPSG:5070")
  sim_sf_ll <- sf::st_transform(sim_sf, 4326)
  ll_coords <- sf::st_coordinates(sim_sf_ll)
  
  df <- data.frame(
    checklist_id = 1:total_n,
    longitude = ll_coords[, 1],
    latitude  = ll_coords[, 2],
    formatted_date = as.Date("2018-06-01") + sample(1:30, total_n, replace=TRUE),
    locality_id = paste0("L", sampled_loc_indices)
  )
  return(df)
}

##########
# 4. MAIN SIMULATION LOOP
##########

all_results_df <- data.frame()
plot_list_sim1 <- list() 

for (sim in 1:n_sims) {
  
  cat(sprintf("\n===========================================\n"))
  cat(sprintf("=== STARTING SIMULATION %d of %d ===\n", sim, n_sims))
  cat(sprintf("===========================================\n"))
  
  # --- A. GENERATE LANDSCAPE ---
  cat("Generating true cellular landscape...\n")
  
  r_cov <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                       xmin=0, xmax=full_grid_dim * res_m, 
                       ymin=0, ymax=full_grid_dim * res_m,
                       vals=scale(rnorm(full_n_cells)))
  names(r_cov) <- "cell_cov1"
  
  # Assign Raster CRS Early
  sf_crs <- "EPSG:5070"
  terra::crs(r_cov) <- sf_crs
  
  fw <- terra::focalMat(r_cov, 3 * res_m, type = "Gauss")
  r_cov <- terra::focal(r_cov, w = fw, fun = sum, na.rm = TRUE)
  terra::values(r_cov) <- as.vector(scale(terra::values(r_cov)))
  
  cellCovs_df <- data.frame(cell_cov1 = terra::values(r_cov, mat=FALSE))
  cellCovs_df$cell_cov1[is.na(cellCovs_df$cell_cov1)] <- 0
  
  X_cell <- model.matrix(~cell_cov1, data = cellCovs_df)
  lambda_j <- exp(X_cell %*% true_betas) * cell_area_km2 
  psi_j <- 1 - exp(-lambda_j)
  Z_j <- rbinom(full_n_cells, 1, psi_j) 
  
  r_abund <- r_cov; terra::values(r_abund) <- lambda_j
  r_occ <- r_cov; terra::values(r_occ) <- Z_j
  
  
  # --- B. SIMULATE OBSERVER DATA ---
  cat("Simulating observation process...\n")
  
  checklists_df <- simulate_ebird_checklists(total_checklists, n_unique_locations, full_grid_dim * res_m)
  checklists_df$obs_cov1 <- rnorm(total_checklists)
  
  # FIX: Project the Lat/Longs back to meters to find the correct Raster Cell
  pts_vect <- terra::vect(checklists_df, geom=c("longitude", "latitude"), crs="EPSG:4326")
  pts_vect_proj <- terra::project(pts_vect, r_cov)
  cell_ids <- terra::cellFromXY(r_cov, terra::crds(pts_vect_proj))
  checklists_df$cell_id <- cell_ids
  
  # Sanitize against edge-case NAs
  checklists_df <- checklists_df[!is.na(checklists_df$cell_id), ]
  
  checklists_df$true_Z <- Z_j[checklists_df$cell_id]
  logit_p <- true_alphas["alpha_int"] + true_alphas["alpha_cov"] * checklists_df$obs_cov1
  checklists_df$p <- plogis(logit_p)
  checklists_df$species_observed <- rbinom(nrow(checklists_df), 1, checklists_df$p * checklists_df$true_Z)
  
  
  # --- C. CLUSTERING & MODELING LOOP ---
  for (method in methods_to_test) {
    cat(sprintf("\n  >>> Applying Method: %s\n", method))
    
    clust_res <- run_clustering_method(method, checklists_df, state_covs = c("cell_cov1"))
    
    if (is.null(clust_res) || is.null(clust_res$data) || nrow(clust_res$data) == 0) {
        cat("      Method failed to cluster or filtered out all data. Skipping.\n")
        next
    }
    
    clust_df <- clust_res$data
    
    cat("      Generating geometries & W matrix...\n")
    site_geoms <- create_site_geometries(clust_df, r_cov, buffer_m, method)
    
    # FIX: Ensure disjoint function knows the input df points are lat/long
    split_res <- disjoint_site_geometries(site_geoms, clust_df, crs_points = 4326)
    
    final_geoms <- split_res$geoms
    final_clust_df <- split_res$data
    
    w_matrix <- generate_overlap_matrix(final_geoms, r_cov)
    
    cat("      Fitting occuPPM...\n")
    umf <- prepare_occuPPM_data(checklists_df, final_clust_df, w_matrix, c("obs_cov1"), cellCovs_df)
    
    fm <- fit_occuPPM_model(
        umf, ~obs_cov1, ~cell_cov1,
        n_reps = n_reps, stable_reps = 5,
        optimizer = selected_optimizer,
        lower = PARAM_LOWER, upper = PARAM_UPPER,
        init_lower = INIT_LOWER, init_upper = INIT_UPPER
    )
    
    if (is.null(fm)) {
        cat("      Model failed to converge.\n")
        est_val <- c(NA, NA, NA, NA)
    } else {
        est_val <- c(coef(fm, 'det'), coef(fm, 'state'))
    }
    
    loop_results <- data.frame(
      Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
      True_Value = c(true_alphas, true_betas),
      Estimated_Value = est_val,
      Method = method,
      sim_id = sim,
      M_sites = nrow(final_geoms)
    )
    all_results_df <- rbind(all_results_df, loop_results)
    
    
    # --- D. PLOTTING FOR SIM 1 ---
    if (sim == 1) {
      cat("      Generating plots for 8x3 grid...\n")
      
      bg_df_cov <- as.data.frame(r_cov, xy=TRUE); names(bg_df_cov)[3] <- "val"
      bg_df_abund <- as.data.frame(r_abund, xy=TRUE); names(bg_df_abund)[3] <- "val"
      bg_df_occ <- as.data.frame(r_occ, xy=TRUE); names(bg_df_occ)[3] <- "val"
      
      geom_sf <- sf::st_as_sf(final_geoms)
      
      # FIX: Plotting points must be projected from Lat/Long to Albers for overlay
      pts_sf <- sf::st_as_sf(final_clust_df, coords=c("longitude", "latitude"), crs=4326)
      pts_sf <- sf::st_transform(pts_sf, sf_crs)
      
      base_theme <- theme_void() + theme(
        plot.title = element_text(size=10, hjust=0.5, face="bold"),
        legend.position = "none"
      )
      
      build_panel <- function(bg_data, fill_scale, title_suffix) {
        ggplot() +
          geom_raster(data = bg_data, aes(x=x, y=y, fill=val)) +
          fill_scale +
          geom_sf(data = geom_sf, fill=NA, color="red", linewidth=0.3) +
          geom_sf(data = pts_sf, color="black", size=0.2, alpha=0.5) +
          coord_sf(expand=FALSE) +
          labs(title = paste(method, "-", title_suffix)) +
          base_theme
      }
      
      p_cov <- build_panel(bg_df_cov, scale_fill_viridis_c(), "Covariate")
      p_abund <- build_panel(bg_df_abund, scale_fill_viridis_c(option="magma"), "Abundance")
      p_occ <- build_panel(bg_df_occ, scale_fill_gradient(low="navyblue", high="yellow"), "Occupancy")
      
      plot_list_sim1[[paste0(method, "_cov")]] <- p_cov
      plot_list_sim1[[paste0(method, "_abund")]] <- p_abund
      plot_list_sim1[[paste0(method, "_occ")]] <- p_occ
    }
    
  } 
} 

##########
# 5. SAVE RESULTS & GENERATE FINAL PLOTS
##########

cat("\n--- Generating Final Outputs ---\n")

write.csv(all_results_df, file.path(output_dir, "clustering_sim_results.csv"), row.names=FALSE)

if (length(plot_list_sim1) > 0) {
  spatial_grid <- patchwork::wrap_plots(plot_list_sim1, ncol = 3, nrow = length(methods_to_test))
  ggsave(file.path(output_dir, "spatial_method_comparison.png"), 
         plot = spatial_grid, width = 12, height = 24, dpi = 300)
}

all_results_df$Error <- all_results_df$Estimated_Value - all_results_df$True_Value
all_results_df$Method <- factor(all_results_df$Method, levels = methods_to_test)

create_err_plot <- function(param_name, title) {
  df_sub <- all_results_df[all_results_df$Parameter == param_name, ]
  ggplot(df_sub, aes(x = Method, y = Error, fill = Method)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, color="red", linetype="dashed", linewidth=1) +
    theme_bw() +
    labs(title = title, x = NULL, y = "Error (Est - True)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

p_err1 <- create_err_plot("beta (state_int)", "State Intercept")
p_err2 <- create_err_plot("beta (state_cov1)", "State Slope")
p_err3 <- create_err_plot("alpha (det_int)", "Observation Intercept")
p_err4 <- create_err_plot("alpha (det_cov1)", "Observation Slope")

combined_error_plot <- (p_err1 | p_err2) / (p_err3 | p_err4)
ggsave(file.path(output_dir, "parameter_error_boxplots.png"), 
       plot = combined_error_plot, width = 12, height = 10, dpi = 300)

cat("--- Simulation Completed Successfully! ---\n")