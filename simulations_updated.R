################################################################
# Updated Simulation Experiments: Clustering Comparison
# Testing MAUP via 4-Scenario Ablation Study
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

# --- The 4 Scenarios ---
scenarios <- c("Idealized", "Ecotone", "Accessibility", "Opportunistic")
n_sims <- 3 # SET TO 100 FOR FULL RUN
n_reps <- 30 # Number of random-start repetitions for model fitting

# --- Landscape Parameters ---
full_grid_dim <- 200
full_n_cells <- full_grid_dim * full_grid_dim
res_m <- 100 # Each cell is 100m x 100m
cell_area_km2 <- (res_m / 1000)^2

# --- Checklist Simulation Parameters ---
total_checklists <- 10000
n_unique_locations <- 5000 

# --- True parameter values ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
# true_betas <- c(beta_int = 1.0, beta_cov = 1.0)
true_betas <- c(beta_int = 4.0, beta_cov = 1.0)

# --- Methods to Test ---
methods_to_test <- c(
  "lat-long", "1to10", "2to10", 
  "0.25-kmSq", "1-kmSq", "4-kmSq", 
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
# 3. HELPER FUNCTIONS
##########

generate_accessibility_surface <- function(r_template, n_roads = 3) {
  ext <- ext(r_template)
  lines_list <- list()
  for(i in 1:n_roads) {
    p1 <- c(runif(1, ext$xmin, ext$xmax), runif(1, ext$ymin, ext$ymax))
    p2 <- c(runif(1, ext$xmin, ext$xmax), runif(1, ext$ymin, ext$ymax))
    lines_list[[i]] <- st_linestring(matrix(c(p1, p2), ncol=2, byrow=TRUE))
  }
  roads_sf <- st_sfc(lines_list, crs = "EPSG:5070")
  roads_vect <- vect(roads_sf)
  
  dist_rast <- distance(r_template, roads_vect)
  scale_factor <- (ext$xmax - ext$xmin) / 5
  prob_rast <- exp(-dist_rast / scale_factor)
  return(prob_rast)
}

simulate_ebird_checklists_updated <- function(total_n, accessibility_raster, scenario) {
  
  # Determine spatial bias based on scenario
  if (scenario %in% c("Idealized", "Ecotone")) {
    probs <- accessibility_raster
    terra::values(probs) <- 1 # Uniform Random
  } else {
    probs <- accessibility_raster # Road-biased
  }
  
  # Explicitly add replace = TRUE to avoid sample size limits
  sampled_pts <- spatSample(probs, size = total_n * 5, method = "weights", xy = TRUE, replace = TRUE)
  
  sampled_pts$loc_idx <- cellFromXY(accessibility_raster, sampled_pts[,c("x", "y")])
  
  # --- GEOMETRY FIX: JITTER TO PREVENT VORONOI CRASH ---
  # GEOS st_voronoi crashes if points are on a perfect mathematical grid.
  # We apply a slight unique jitter (up to 10m) to each location to break grid degeneracy,
  # but ensure all checklists at the same 'loc_idx' still get the EXACT same jitter.
  uniq_locs <- unique(sampled_pts$loc_idx)
  jit_x <- runif(length(uniq_locs), -10, 10)
  jit_y <- runif(length(uniq_locs), -10, 10)
  match_idx <- match(sampled_pts$loc_idx, uniq_locs)
  
  sampled_pts$x <- sampled_pts$x + jit_x[match_idx]
  sampled_pts$y <- sampled_pts$y + jit_y[match_idx]
  # -----------------------------------------------------
  
  df_temp <- sampled_pts %>%
    group_by(loc_idx) %>%
    mutate(visit_num = row_number()) %>%
    filter(visit_num <= 15) %>% 
    ungroup() %>%
    slice_head(n = total_n)
  
  # Inject Positional Blurring for the Opportunistic scenario
  if (scenario == "Opportunistic") {
    # Observers drop the pin at (x,y), but walked up to 500m away
    angle <- runif(nrow(df_temp), 0, 2 * pi)
    radius <- sqrt(runif(nrow(df_temp), 0, 1)) * 500
    df_temp$true_x <- df_temp$x + radius * cos(angle)
    df_temp$true_y <- df_temp$y + radius * sin(angle)
  } else {
    # Exact Locations
    df_temp$true_x <- df_temp$x
    df_temp$true_y <- df_temp$y
  }
  
  # Convert *reported* coords to Lat/Long (This is what the clustering algorithms see)
  pts_sf <- st_as_sf(df_temp, coords = c("x", "y"), crs = "EPSG:5070")
  pts_sf_ll <- st_transform(pts_sf, 4326)
  ll_coords <- st_coordinates(pts_sf_ll)
  
  df <- data.frame(
    checklist_id = 1:nrow(df_temp),
    longitude = ll_coords[, 1], # Reported Lat/Long
    latitude  = ll_coords[, 2],
    true_x = df_temp$true_x,    # True Albers X
    true_y = df_temp$true_y,    # True Albers Y
    locality_id = paste0("L", df_temp$loc_idx),
    observation_date = as.Date("2026-01-01") + sample(0:30, nrow(df_temp), replace = TRUE)
  )
  
  # Keep the formatted_date fix for auk::filter_repeat_visits
  df$formatted_date <- df$observation_date 
  
  return(df)
}

fit_occuPPM_model_sim <- function(umf, state_formula, obs_formula, n_reps = 30, 
                                  stable_reps = 10, optimizer = "nlminb", 
                                  lower = -Inf, upper = Inf,
                                  init_lower = -Inf, init_upper = Inf) {
  
  occuPPM_formula <- as.formula(paste(
    paste(deparse(obs_formula), collapse = ""), 
    paste(deparse(state_formula), collapse = "")
  ))
  
  n_obs_pars <- length(all.vars(obs_formula)) + 1 
  n_state_pars <- length(all.vars(state_formula)) + 1
  n_params <- n_obs_pars + n_state_pars
  
  best_fm <- NULL
  min_nll <- Inf
  fit_successful <- FALSE
  stable_count <- 0
  tolerance <- 0.01 
  
  for (rep in 1:n_reps) {
    rand_starts <- runif(n_params, min = init_lower, max = init_upper)
    
    fm_rep <- try(unmarked::occuPPM(
      formula = occuPPM_formula, data = umf, starts = rand_starts,
      se = FALSE, method = optimizer, lower = lower, upper = upper
    ), silent = TRUE)
    
    if (!inherits(fm_rep, "try-error")) {
      current_nll <- fm_rep@negLogLike
      if (is.finite(current_nll)) {
        if (current_nll < min_nll) {
          if (abs(min_nll - current_nll) < tolerance) { stable_count <- stable_count + 1 } else { stable_count <- 0 }
          min_nll <- current_nll
          best_fm <- fm_rep
          fit_successful <- TRUE
        } else if (abs(current_nll - min_nll) < tolerance) {
          stable_count <- stable_count + 1
        }
        if (stable_count >= stable_reps) break 
      }
    }
  }
  if (!fit_successful) return(NULL)
  return(best_fm)
}

##########
# 4. MAIN SIMULATION LOOP
##########

all_results_df <- data.frame()
plot_collection <- list()
landscape_plots <- list() # For the 4-panel landscape plot

for (scen in scenarios) {
  plot_collection[[scen]] <- list()
  
  for (sim in 1:n_sims) {
    
    cat(sprintf("\n===========================================\n"))
    cat(sprintf("=== SCENARIO: %s | SIM %d of %d ===\n", scen, sim, n_sims))
    cat(sprintf("===========================================\n"))
    
    # --- A. GENERATE LANDSCAPE ---
    cat("Generating true cellular landscape...\n")
    
    sf_crs <- "EPSG:5070"
    
    # 1. Base Noise
    r_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                           xmin=0, xmax=full_grid_dim * res_m, 
                           ymin=0, ymax=full_grid_dim * res_m,
                           vals=rnorm(full_n_cells))
    terra::crs(r_noise) <- sf_crs
    
    # 2. Smooth Base
    fw <- terra::focalMat(r_noise, 15 * res_m, type = "Gauss")
    r_smooth <- terra::focal(r_noise, w = fw, fun = mean, na.rm = TRUE)
    r_smooth <- terra::scale(r_smooth)
    
    access_rast <- generate_accessibility_surface(r_noise)
    
    # 3. Apply Thresholding for Ecotone Scenarios
    if (scen == "Idealized") {
      r_cov <- r_smooth
    } else {
      # Classify into 3 distinct habitat patches
      r_habitat <- terra::classify(r_smooth, 
                                   rcl=matrix(c(-Inf, -0.5, -1.5, 
                                                -0.5, 0.5,  0.0, 
                                                 0.5, Inf,  1.5), ncol=3, byrow=TRUE))
      # Add small local variance back in
      r_micro <- r_noise
      terra::values(r_micro) <- rnorm(full_n_cells, sd=0.3)
      r_cov <- r_habitat + r_micro
    }
    
    names(r_cov) <- "cell_cov1"
    
    cellCovs_df <- data.frame(cell_cov1 = terra::values(r_cov, mat=FALSE))
    cellCovs_df$cell_cov1[is.na(cellCovs_df$cell_cov1)] <- 0
    
    X_cell <- model.matrix(~cell_cov1, data = cellCovs_df)
    lambda_j <- exp(X_cell %*% true_betas) * cell_area_km2 
    psi_j <- 1 - exp(-lambda_j)
    Z_j <- rbinom(full_n_cells, 1, psi_j) 
    
    # --- B. SIMULATE OBSERVER DATA ---
    cat("Simulating observation process...\n")
    
    checklists_df <- simulate_ebird_checklists_updated(5000, access_rast, scenario = scen)
    checklists_df$obs_cov1 <- rnorm(nrow(checklists_df))
    
    # IMPORTANT: Extract true biology based on TRUE location, not reported pin
    pts_true_vect <- terra::vect(as.matrix(checklists_df[, c("true_x", "true_y")]), type="points", crs=sf_crs)
    cell_ids_true <- terra::cellFromXY(r_cov, terra::crds(pts_true_vect))
    checklists_df$cell_id <- cell_ids_true
    checklists_df <- checklists_df[!is.na(checklists_df$cell_id), ]
    
    checklists_df$true_Z <- Z_j[checklists_df$cell_id]
    logit_p <- true_alphas["alpha_int"] + true_alphas["alpha_cov"] * checklists_df$obs_cov1
    checklists_df$p <- plogis(logit_p)
    checklists_df$species_observed <- rbinom(nrow(checklists_df), 1, checklists_df$p * checklists_df$true_Z)
    
    # Store true covariate for later comparison
    checklists_df$cell_cov1 <- cellCovs_df$cell_cov1[checklists_df$cell_id]
    
    # --- Generate Landscape Plot for Sim 1 ---
    if (sim == 1) {
      bg_df <- as.data.frame(r_cov, xy=TRUE)
      p_land <- ggplot() +
        geom_raster(data = bg_df, aes(x=x, y=y, fill=cell_cov1)) +
        scale_fill_viridis_c() +
        geom_point(data = checklists_df, aes(x=true_x, y=true_y), size=0.1, alpha=0.3, color="red") +
        coord_fixed(expand=FALSE) + theme_void() +
        theme(legend.position="none", plot.title=element_text(hjust=0.5, face="bold")) +
        labs(title = scen)
      
      landscape_plots[[scen]] <- p_land
    }
    
    # Count visits per unique *reported* location
    checklists_df <- checklists_df %>%
      dplyr::group_by(locality_id) %>%
      dplyr::mutate(n_visits = dplyr::n()) %>%
      dplyr::ungroup()
    
    checklists_1to10_df <- checklists_df %>% dplyr::filter(n_visits >= 1 & n_visits <= 10)
    checklists_2to10_df <- checklists_df %>% dplyr::filter(n_visits >= 2 & n_visits <= 10)
    
    # --- C. CLUSTERING & MODELING LOOP ---
    for (method in methods_to_test) {
      cat(sprintf("\n  >>> Applying Method: %s\n", method))
      
      if (method == "lat-long") { current_df <- checklists_df }
      else if (method == "1to10") { current_df <- checklists_1to10_df }
      else { current_df <- checklists_2to10_df }
      
      clust_res <- run_clustering_method(method, current_df, state_covs = c("cell_cov1"))
      if (is.null(clust_res) || is.null(clust_res$data) || nrow(clust_res$data) == 0) next
      
      clust_df <- clust_res$data
      
      site_geoms <- create_site_geometries(clust_df, r_cov, buffer_m, method)
      split_res <- disjoint_site_geometries(site_geoms, clust_df, crs_points = 4326)
      final_geoms <- split_res$geoms
      final_clust_df <- split_res$data
      
      w_matrix <- generate_overlap_matrix(final_geoms, r_cov)
      umf <- prepare_occuPPM_data(current_df, final_clust_df, w_matrix, c("obs_cov1"), cellCovs_df)
      
      fm <- fit_occuPPM_model_sim(
        umf, state_formula = ~cell_cov1, obs_formula = ~obs_cov1,
        n_reps = n_reps, stable_reps = 5, optimizer = selected_optimizer,
        lower = PARAM_LOWER, upper = PARAM_UPPER, init_lower = INIT_LOWER, init_upper = INIT_UPPER
      )
      
      est_val <- if (is.null(fm)) c(NA, NA, NA, NA) else c(coef(fm, 'det'), coef(fm, 'state'))
      
      loop_results <- data.frame(
        Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
        True_Value = c(true_alphas, true_betas),
        Estimated_Value = est_val,
        Method = method, Scenario = scen, sim_id = sim, M_sites = nrow(final_geoms)
      )
      all_results_df <- rbind(all_results_df, loop_results)
      
      if (sim == 1) {
        geom_sf <- sf::st_as_sf(final_geoms)
        pts_sf <- sf::st_as_sf(final_clust_df, coords=c("longitude", "latitude"), crs=4326)
        pts_sf <- sf::st_transform(pts_sf, sf_crs)
        
        p_method <- ggplot() +
          geom_raster(data = as.data.frame(r_cov, xy=TRUE), aes(x=x, y=y, fill=cell_cov1), alpha=0.3) +
          geom_sf(data = geom_sf, fill=NA, color="red", linewidth=0.3) +
          geom_sf(data = pts_sf, color="black", size=0.1, alpha=0.3) +
          coord_sf(expand=FALSE) + labs(title = method) +
          theme_void() + theme(legend.position="none", plot.title=element_text(size=8))
        plot_collection[[scen]][[method]] <- p_method
      }
    } 
  }
} 

##########
# 5. SAVE RESULTS & GENERATE FINAL PLOTS
##########

cat("\n--- Generating Final Outputs ---\n")
write.csv(all_results_df, file.path(output_dir, "clustering_sim_results.csv"), row.names=FALSE)

# --- A. Save Landscape Overview Plot (4 Panels) ---
landscape_grid <- patchwork::wrap_plots(landscape_plots, ncol = 2) + 
  plot_annotation(title = "Simulation Scenarios (True Covariate + Sampling Locations)")
ggsave(file.path(output_dir, "landscape_scenarios.png"), plot = landscape_grid, width = 10, height = 10, dpi = 300)

# --- B. Parameter Error Boxplots (Facetted by 4 Scenarios) ---
all_results_df$Error <- all_results_df$Estimated_Value - all_results_df$True_Value
all_results_df$Method <- factor(all_results_df$Method, levels = methods_to_test)
# Force factor order so facets display logically
all_results_df$Scenario <- factor(all_results_df$Scenario, levels = scenarios)

create_err_plot <- function(param_name, title) {
  df_sub <- all_results_df[all_results_df$Parameter == param_name, ]
  
  ggplot(df_sub, aes(x = Method, y = Error, fill = Method)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, color="red", linetype="dashed", linewidth=1) +
    facet_grid(~ Scenario) + 
    theme_bw() +
    labs(title = title, x = NULL, y = "Error (Est - True)") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "bottom",
      strip.background = element_rect(fill="lightgrey")
    )
}

p_err1 <- create_err_plot("beta (state_int)", "State Intercept")
p_err2 <- create_err_plot("beta (state_cov1)", "State Slope")
p_err3 <- create_err_plot("alpha (det_int)", "Observation Intercept")
p_err4 <- create_err_plot("alpha (det_cov1)", "Observation Slope")


combined_error_plot <- (p_err1 | p_err2) / (p_err3 | p_err4) + 
  plot_layout(guides = "collect")

ggsave(file.path(output_dir, "parameter_error_boxplots.png"), plot = combined_error_plot, width = 18, height = 14, dpi = 300)

cat("--- Simulation Completed Successfully! ---\n")