#################################################################
# Simulation Experiments using clustGeo Sites
# Decoupled Data Generating Process: 
# Sites defined organically via ClustGeo + Voronoi, uniformly sampled (J=3)
#################################################################

###
# 1. SETUP
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
library(scales) 
library(sf)
library(dplyr)
library(tidyr)
library(ClustGeo)

# Source required helpers
source(file.path("R", "utils.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Simulation repetitions ---
n_sims <- 3
n_reps <- 30 

# --- Landscape Configuration ---
full_grid_dim <- 200 
full_n_cells <- full_grid_dim * full_grid_dim # 40000

# --- ClustGeo Sampling Frame Configurations ---
dummy_N_value <- 4800          # Anchor points used to draw geometries (NOT the final visits)
M_total_sites <- 1600          # Total number of contiguous territories to form
visits_per_site <- 3           # J=3 virtual observers per site (eliminates area-visit confounding)
clustgeo_alpha <- 0.5          # Blends environment and spatial distance to define territory shapes

# --- True parameter values (Restored to optimal simulation biology) ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
true_betas <- c(beta_int = -5.0, beta_cov = 1.0)

# --- Model settings ---
selected_optimizer <- "nlminb"
PARAM_LOWER <- -20
PARAM_UPPER <- 20

# --- Ablation Study Parameters ---
M_values_to_test <- c(100, 200, 400, 800, 1600)

# --- Spatial Autocorrelation (SAC) Settings ---
sac_levels <- c("Low", "Medium", "High") 
sac_sigmas <- c(Low = 0, Medium = 3, High = 6) 

# --- Skew Patterns ---
skew <- "Centers"
n_centers <- 1
centers_scale <- 5
decay_scale <- 30^2

output_dir <- file.path("output", "simulation_experiments", "clustgeo")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- Simulation Starting ---\n")
cat(sprintf("Running %d full simulations per SAC level.\n", n_sims))
cat(sprintf("Dummy Anchors: %d | Clusters: %d | Alpha: %.1f\n", dummy_N_value, M_total_sites, clustgeo_alpha))
cat(sprintf("SAC Levels: %s\n", paste(sac_levels, collapse=", ")))

##########
# 3. Local Spatial Helpers (Coordinate-Agnostic)
##########

# Generates W Matrix without complaining about meters or lat/long CRS mismatches
local_generate_overlap_matrix <- function(site_geoms_sf, reference_raster) {
  site_geoms_sf$site <- as.character(site_geoms_sf$site)
  site_vect <- terra::vect(site_geoms_sf)
  
  overlap_df <- terra::extract(reference_raster[[1]], site_vect, cells = TRUE, exact = TRUE, ID = TRUE)
  
  val_col_name <- names(reference_raster)[1]
  overlap_df <- overlap_df[!is.na(overlap_df[[val_col_name]]), ]
  
  if (!("fraction" %in% names(overlap_df))) {
    overlap_df$fraction <- 1.0
  }
  overlap_df$w_area <- overlap_df$fraction 
  
  n_sites <- nrow(site_geoms_sf)
  n_cells <- terra::ncell(reference_raster)
  
  w <- Matrix::sparseMatrix(
    i = overlap_df$ID,     
    j = overlap_df$cell,   
    x = overlap_df$w_area,
    dims = c(n_sites, n_cells),
    dimnames = list(site_geoms_sf$site, NULL) 
  )
  return(w)
}

##########
# 4. Pre-generate Fixed Seeds
##########

full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1

cat("Pre-generating fixed COVARIATE skew seeds...\n")
cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(n_centers, 0, full_grid_dim), 
    y = runif(n_centers, 0, full_grid_dim)
  )
}

results_storage <- list()

##########
# 5. OUTER LOOPS: SAC -> SIM
##########

for (sac_level in sac_levels) {
  
  current_sigma <- sac_sigmas[[sac_level]]
  
  for (sim in 1:n_sims) {
    
    cat(sprintf("\n=== SAC: %s | Sim %d of %d ===\n", sac_level, sim, n_sims))
    
    # --- 5a. Generate Base Spatial Noise (SAC) ---
    r_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                           xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim,
                           vals=rnorm(full_n_cells))
    
    if (current_sigma > 0) {
      fw <- terra::focalMat(r_noise, current_sigma, type = "Gauss")
      r_smooth <- terra::focal(r_noise, w = fw, fun = sum, na.rm = TRUE)
    } else {
      r_smooth <- r_noise
    }
    terra::values(r_smooth) <- as.vector(scale(terra::values(r_smooth)))
    
    # --- 5b. Apply Skew / Trend ---
    seeds <- cov_center_seeds[[sim]] 
    
    r_centers <- terra::rast(r_smooth)
    terra::values(r_centers) <- 0
    
    rows <- terra::init(r_smooth, "y")
    cols <- terra::init(r_smooth, "x")
    
    for(k in seq_along(seeds$x)) {
      d2 <- (cols - seeds$x[k])^2 + (rows - seeds$y[k])^2
      r_centers <- r_centers + exp(-d2 / (2 * decay_scale))
    }
    
    r_trend <- r_smooth + (r_centers * centers_scale)
    terra::values(r_trend) <- as.vector(scale(terra::values(r_trend)))
    r_final <- r_trend
    names(r_final) <- "cell_cov1"
    
    full_cellCovs <- data.frame(cell_cov1 = terra::values(r_final, mat=FALSE))
    if(any(is.na(full_cellCovs$cell_cov1))) full_cellCovs$cell_cov1[is.na(full_cellCovs$cell_cov1)] <- 0
    
    full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
    full_lambda_j <- exp(full_X_cell %*% true_betas)
    
    # Generate pixel-level true occupancy based on abundance (for plotting)
    full_pixel_psi <- 1 - exp(-full_lambda_j)
    full_pixel_Z <- rbinom(full_n_cells, 1, full_pixel_psi)
    
    # --- 5c. PHASE 1: Create Dummy Anchor Points for Geometry ---
    dummy_x <- runif(dummy_N_value, 0, full_grid_dim)
    dummy_y <- runif(dummy_N_value, 0, full_grid_dim)
    dummy_df <- data.frame(id = 1:dummy_N_value, x = dummy_x, y = dummy_y)
    
    dummy_vect <- terra::vect(dummy_df, geom=c("x", "y"))
    dummy_df$cell_cov1 <- terra::extract(r_final, dummy_vect)$cell_cov1
    
    # ClustGeo grouping of the dummy anchors
    env_dist <- dist(scale(dummy_df$cell_cov1))
    geo_dist <- dist(scale(dummy_df[, c("x", "y")]))
    tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = clustgeo_alpha)
    dummy_df$site <- cutree(tree, M_total_sites)
    
    dummy_sf <- sf::st_as_sf(dummy_df, coords=c("x", "y"), remove=FALSE)
    
    # --- 5d. PHASE 2: Voronoi Tessellation & Dissolve (Creates Puzzle Pieces) ---
    landscape_boundary <- sf::st_polygon(list(matrix(c(
      0, 0, 
      full_grid_dim, 0, 
      full_grid_dim, full_grid_dim, 
      0, full_grid_dim, 
      0, 0
    ), ncol = 2, byrow = TRUE))) %>%
      sf::st_sfc() %>%
      sf::st_sf(geometry = .)
    
    bbox_polygon <- sf::st_as_sfc(sf::st_bbox(landscape_boundary))
    
    # Generate gapless grid of Voronoi tiles
    voronoi_tiles <- sf::st_voronoi(sf::st_union(dummy_sf), envelope = bbox_polygon, dTolerance = 0) %>%
      sf::st_collection_extract(type = "POLYGON") %>%
      sf::st_sf()
    
    # Merge tiles based on ClustGeo ID to form large, non-overlapping territories
    voronoi_w_id <- sf::st_join(voronoi_tiles, dummy_sf, join = sf::st_contains)
    
    final_site_geoms <- voronoi_w_id %>%
      dplyr::group_by(site) %>%
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
      sf::st_make_valid()
    
    final_site_geoms <- suppressWarnings(sf::st_intersection(final_site_geoms, landscape_boundary))
    
    w_matrix <- local_generate_overlap_matrix(final_site_geoms, r_final)
    
    # --- 5e. PHASE 3: True Biology and Standardized J=3 Observers ---
    lambda_tilde_i <- as.numeric(w_matrix %*% full_lambda_j)
    psi_i <- 1 - exp(-lambda_tilde_i)
    Z_i <- rbinom(length(psi_i), 1, psi_i)
    site_Z <- setNames(Z_i, rownames(w_matrix))
    
    # Decoupling: Drop exactly 3 uniform observation points inside every site
    pts_list <- lapply(1:nrow(final_site_geoms), function(i) {
      site_geom <- final_site_geoms[i, ]
      pts <- suppressMessages(sf::st_sample(site_geom, size = visits_per_site, type = "random", exact = TRUE))
      
      # Safety fallback for exceptionally tiny polygons
      if (length(pts) < visits_per_site) {
        centroid <- sf::st_centroid(site_geom)
        pts <- c(pts, rep(centroid$geometry, visits_per_site - length(pts)))
      } else if (length(pts) > visits_per_site) {
        pts <- pts[1:visits_per_site]
      }
      sf::st_sf(site = as.character(site_geom$site), geometry = pts)
    })
    
    final_points <- do.call(rbind, pts_list)
    final_points$obs_cov1 <- rnorm(nrow(final_points))
    final_points$Z_i <- site_Z[final_points$site]
    
    logit_p <- true_alphas["alpha_int"] + true_alphas["alpha_cov"] * final_points$obs_cov1
    p <- plogis(logit_p)
    final_points$y <- rbinom(nrow(final_points), 1, p * final_points$Z_i)
    
    # Initialize plotting list for this sim
    if (sim == 1) {
      plots_cov   <- list()
      plots_abund <- list()
      plots_occ   <- list()
    }
    
    # --- 5f. Loop Over M (Sample Size) ---
    total_sites <- nrow(w_matrix)
    valid_M_values <- M_values_to_test[M_values_to_test <= total_sites]
    
    for (M_i in valid_M_values) {
      
      # Uniformly select sites
      sel_sites <- sample(rownames(w_matrix), M_i, replace = FALSE)
      sel_sites <- sort(sel_sites)
      
      # Subset Data
      sub_w <- w_matrix[sel_sites, , drop=FALSE]
      sub_pts <- final_points %>% filter(site %in% sel_sites)
      
      # Prepare Format for occuPPM
      occu_df <- sf::st_drop_geometry(sub_pts) %>%
        group_by(site) %>%
        mutate(visit = row_number()) %>%
        ungroup()
      
      # PERFECT ROW ALIGNMENT: Match pivoted data to the exact row order of sub_w
      y_df <- occu_df %>% dplyr::select(site, visit, y) %>%
        tidyr::pivot_wider(names_from=visit, values_from=y)
      y_df <- y_df[match(rownames(sub_w), y_df$site), ] 
      y_wide <- y_df %>% dplyr::select(-site) %>% as.matrix()
      
      obs_cov1_df <- occu_df %>% dplyr::select(site, visit, obs_cov1) %>%
        tidyr::pivot_wider(names_from=visit, values_from=obs_cov1)
      obs_cov1_df <- obs_cov1_df[match(rownames(sub_w), obs_cov1_df$site), ] 
      obs_cov1_wide <- obs_cov1_df %>% dplyr::select(-site) %>% as.matrix()
      
      # Fit Model
      umf <- unmarkedFrameOccuPPM(
        y = y_wide,
        obsCovs = list(obs_cov1 = obs_cov1_wide),
        cellCovs = full_cellCovs, 
        w = sub_w 
      )
      
      best_fm <- NULL
      min_nll <- Inf
      n_params <- length(true_alphas) + length(true_betas)
      
      for (rep in 1:n_reps) {
        rand_starts <- runif(n_params, -6, 2) 
        fm_rep <- try(suppressWarnings(occuPPM(
          formula = ~obs_cov1 ~ cell_cov1,
          data = umf,
          starts = rand_starts,
          se = TRUE,
          method = selected_optimizer,
          lower = PARAM_LOWER, 
          upper = PARAM_UPPER 
        )), silent = TRUE)
        
        if (inherits(fm_rep, "try-error")) next
        
        current_nll <- fm_rep@negLogLike
        if (current_nll < min_nll) {
          min_nll <- current_nll
          best_fm <- fm_rep
        }
      } 
      
      fm <- best_fm
      est_val <- if(is.null(fm)) c(NA,NA,NA,NA) else c(coef(fm, 'det'), coef(fm, 'state'))
      
      # Record Result
      loop_results <- data.frame(
        Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
        True_Value = c(true_alphas, true_betas),
        Estimated_Value = est_val,
        M = M_i, 
        sim_id = sim,
        SAC_Level = sac_level
      )
      
      results_storage[[length(results_storage) + 1]] <- loop_results
      
      # --- PLOTTING LOGIC (Inside M loop) ---
      if (sim == 1) {
        
        cell_df <- data.frame(
          x = full_cell_col,
          y = full_cell_row,
          covariate = full_cellCovs$cell_cov1,
          pixel_abundance = full_lambda_j,
          pixel_occupancy = as.factor(full_pixel_Z)
        )
        
        sel_geoms <- final_site_geoms %>% filter(site %in% sel_sites)
        
        tight_theme <- theme_minimal() + 
          theme(
            axis.title = element_blank(),
            plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt"),
            legend.position = "bottom",
            legend.direction = "horizontal"
          )
        
        # Col 1: Covariate + Red Borders + Black Points
        p_cov <- ggplot(cell_df, aes(x=x, y=y, fill=covariate)) +
          geom_raster() +
          scale_fill_viridis_c() +
          geom_sf(data=sel_geoms, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          geom_point(data=sf::st_drop_geometry(sub_pts), aes(x=x, y=y), color="black", size=0.5, inherit.aes=FALSE) +
          coord_sf(expand=FALSE, datum=NA) + 
          labs(title=sprintf("Covariate (M=%d)", M_i), fill="Covariate") +
          tight_theme
        
        # Col 2: Abundance + Red Borders
        p_abund <- ggplot(cell_df, aes(x=x, y=y, fill=pixel_abundance)) +
          geom_raster() +
          scale_fill_viridis_c(option = "magma") +
          geom_sf(data=sel_geoms, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          coord_sf(expand=FALSE, datum=NA) + 
          labs(title=sprintf("Abundance (M=%d)", M_i), fill="Abundance") +
          tight_theme
        
        # Col 3: Occupancy + Red Borders
        p_occ <- ggplot(cell_df, aes(x=x, y=y, fill=pixel_occupancy)) +
          geom_raster() +
          scale_fill_manual(values=c("0"="navyblue", "1"="yellow")) +
          geom_sf(data=sel_geoms, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          coord_sf(expand=FALSE, datum=NA) + 
          labs(title=sprintf("Occupancy (M=%d)", M_i), fill="Occupancy") +
          tight_theme
          
        plots_cov[[length(plots_cov)+1]]     <- p_cov
        plots_abund[[length(plots_abund)+1]] <- p_abund
        plots_occ[[length(plots_occ)+1]]     <- p_occ
      }
      
    } # End M Loop
    
    # --- SAVE PLOTS ---
    if (sim == 1) {
      cat(sprintf("\nSaving plots for SAC=%s...\n", sac_level))
      
      col_cov <- patchwork::wrap_plots(plots_cov, ncol = 1) 
      col_abund <- patchwork::wrap_plots(plots_abund, ncol = 1) 
      col_occ <- patchwork::wrap_plots(plots_occ, ncol = 1)
      
      final_comb_plot <- (col_cov | col_abund | col_occ) + patchwork::plot_layout(guides = "collect")
      
      fname <- sprintf("plot_SAC=%s.png", sac_level)
      
      suppressMessages(suppressWarnings({
        ggsave(file.path(output_dir, fname), plot=final_comb_plot, dpi=300, width=11, height=20)
      }))
    }
    
    gc()
  } # End Sim Loop
} # End SAC Loop

##########
# 6. Save Results & Generate Error Plots 
##########

cat("\n--- All Simulations Complete. Saving Results... ---\n")

results_df <- do.call(rbind, results_storage)

if (!is.null(results_df) && nrow(results_df) > 0) {
    
    results_df <- results_df[!sapply(results_df$Parameter, is.null), ]
    
    # Save CSV
    write.csv(results_df, file.path(output_dir, "params_clustgeo.csv"), row.names = FALSE)
    
    # Boxplots
    strat_df <- results_df
    strat_df$Error <- strat_df$True_Value - strat_df$Estimated_Value
    strat_df$M_factor <- as.factor(strat_df$M)
    strat_df$SAC_Level <- factor(strat_df$SAC_Level, levels = c("Low", "Medium", "High"))
    
    create_error_plot <- function(param_name, title) {
      ggplot(strat_df[strat_df$Parameter == param_name, ], 
             aes(x = M_factor, y = Error, fill = SAC_Level)) +
        geom_boxplot(outlier.size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_fill_manual(values = c("Low" = "yellow", "Medium" = "orange", "High" = "red")) +
        labs(title = title, x = "M (Sites)", y = "Error (True - Estimate)", fill = "Spatial Autocorrelation") +
        theme_bw() + 
        theme(legend.position = "bottom", legend.direction = "horizontal") 
    }
    
    p1 <- create_error_plot("beta (state_int)", "State Intercept")
    p2 <- create_error_plot("beta (state_cov1)", "State Slope")
    p3 <- create_error_plot("alpha (det_int)", "Observation Intercept")
    p4 <- create_error_plot("alpha (det_cov1)", "Observation Slope")
    
    # Combine natively without the problematic '&' operator
    combined_error_plot <- (p1 | p2) / (p3 | p4) + patchwork::plot_layout(guides = "collect")
    
    suppressMessages(suppressWarnings({
        ggsave(file.path(output_dir, "error_boxplots.png"), plot = combined_error_plot, dpi = 300, width = 10, height = 10)
    }))
    
    cat("\nSaved results and error boxplots.")
}

cat("\n--- Script Finished Successfully ---\n")