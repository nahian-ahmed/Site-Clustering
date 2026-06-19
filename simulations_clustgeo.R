#################################################################
# Simulation Experiments using clustGeo Sites
# Spatial Clustering with uniform sampling and cell-level visualizations
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

# --- ClustGeo & Point Generation Configurations ---
max_N_value <- 4800            # Total observation points generated 
kappa_for_clustgeo <- 33.3     # 33.3% of 4800 = ~1600 sites (Average 3 visits per site)

buffer_cells <- 2              # Buffer in arbitrary grid units (not meters)

# --- True parameter values ---
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
cat(sprintf("Points generated: %d | ClustGeo Kappa: %.1f\n", max_N_value, kappa_for_clustgeo))
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
  
  # Because each cell is 1x1, area is literally just the fraction
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

# Processes disjoint shapes safely in our un-projected spatial frame
local_disjoint_site_geometries <- function(site_geoms_sf, points_sf) {
  sites_split <- site_geoms_sf %>%
    sf::st_make_valid() %>%
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>% 
    sf::st_cast("POLYGON", warn = FALSE)
  
  if (nrow(sites_split) == nrow(site_geoms_sf)) {
    return(list(geoms = site_geoms_sf, data = points_sf))
  }
  
  sites_split <- sites_split %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(sub_id = dplyr::row_number(), new_site_id = paste0(site, "_", sub_id)) %>%
    dplyr::ungroup()
  
  join_res <- sf::st_join(points_sf, sites_split["new_site_id"], join = sf::st_intersects, left = TRUE)
  join_res <- join_res[!duplicated(join_res$checklist_id), ]
  
  points_sf$site <- ifelse(!is.na(join_res$new_site_id), join_res$new_site_id, as.character(points_sf$site))
  
  occupied_new_ids <- unique(na.omit(points_sf$site))
  
  sites_final <- sites_split %>%
    dplyr::filter(new_site_id %in% occupied_new_ids) %>%
    dplyr::select(-site, -sub_id) %>%
    dplyr::rename(site = new_site_id) %>%
    dplyr::select(site, geometry)
  
  return(list(geoms = sites_final, data = points_sf))
}

# Generates Non-Overlapping shapes scaled for a 200x200 grid
local_voronoi_clipped_buffers <- function(points_sf, buffer_dist) {
  bbox_polygon <- sf::st_as_sfc(sf::st_bbox(points_sf) + buffer_dist * 3)
  
  voronoi_tiles <- sf::st_voronoi(sf::st_union(points_sf), envelope = bbox_polygon, dTolerance = 0) %>%
    sf::st_collection_extract(type = "POLYGON") %>%
    sf::st_sf()
  
  voronoi_w_id <- sf::st_join(voronoi_tiles, points_sf, join = sf::st_contains)
  
  site_territories <- voronoi_w_id %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    sf::st_make_valid()
  
  site_buffers <- points_sf %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(geometry = sf::st_convex_hull(sf::st_union(geometry)), .groups = "drop") %>%
    sf::st_buffer(dist = buffer_dist) %>%
    sf::st_make_valid()
  
  site_territories <- dplyr::rename(site_territories, site_t = site)
  site_buffers <- dplyr::rename(site_buffers, site_b = site)
  
  sf::st_agr(site_territories) <- "constant"
  sf::st_agr(site_buffers) <- "constant"
  
  intersections <- suppressWarnings(sf::st_intersection(site_territories, site_buffers))
  
  final_geoms <- intersections %>%
    dplyr::filter(site_t == site_b) %>%
    dplyr::rename(site = site_t) %>%
    dplyr::select(site, geometry) %>%
    sf::st_make_valid() 
  
  # THE FIX: Scale dTolerance down for a 200x200 coordinate grid
  final_geoms <- sf::st_simplify(final_geoms, dTolerance = 0.1, preserveTopology = TRUE)

  geom_type <- sf::st_geometry_type(final_geoms, by_geometry = FALSE)
  if (inherits(geom_type, "GEOMETRYCOLLECTION") || any(geom_type == "GEOMETRYCOLLECTION")) {
      final_geoms <- sf::st_collection_extract(final_geoms, "POLYGON")
  }

  return(final_geoms)
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
    
    # --- 5c. Generate Observation Points & ClustGeo ---
    pts_x <- runif(max_N_value, 0, full_grid_dim)
    pts_y <- runif(max_N_value, 0, full_grid_dim)
    pts_df <- data.frame(checklist_id = 1:max_N_value, x = pts_x, y = pts_y)
    
    pts_vect <- terra::vect(pts_df, geom=c("x", "y"))
    pts_df$cell_cov1 <- terra::extract(r_final, pts_vect)$cell_cov1
    
    # Clustering
    env_dist <- dist(scale(pts_df$cell_cov1))
    geo_dist <- dist(scale(pts_df[, c("x", "y")]))
    
    tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = 1.0)
    num_clusters <- max(2, round(max_N_value * (kappa_for_clustgeo / 100.0)))
    pts_df$site <- cutree(tree, num_clusters)
    
    # --- 5d. Polygonize and Generate W Matrix ---
    points_sf <- sf::st_as_sf(pts_df, coords=c("x", "y"), remove=FALSE)
    
    # Generate geometries using our local scale-adjusted function
    site_geoms_sf <- local_voronoi_clipped_buffers(points_sf, buffer_dist = buffer_cells)
    
    # NEW: Create a hard boundary for the 200x200 landscape
    landscape_boundary <- sf::st_polygon(list(matrix(c(
      0, 0, 
      full_grid_dim, 0, 
      full_grid_dim, full_grid_dim, 
      0, full_grid_dim, 
      0, 0
    ), ncol = 2, byrow = TRUE))) %>%
      sf::st_sfc() %>%
      sf::st_sf(geometry = .)
    
    # FIX: Repair microscopic topological anomalies before intersecting
    site_geoms_sf <- sf::st_make_valid(site_geoms_sf)
    
    # NEW: Clip all site geometries (silencing the spatial constant warning)
    site_geoms_sf <- suppressWarnings(sf::st_intersection(site_geoms_sf, landscape_boundary))
    
    # Splitting disjoint shapes into independent polygons
    disjoint_res <- local_disjoint_site_geometries(site_geoms_sf, points_sf)
    final_geoms <- disjoint_res$geoms
    final_points <- disjoint_res$data
    
    # W Matrix
    w_matrix <- local_generate_overlap_matrix(final_geoms, r_final)
    
    # --- 5e. Simulate True Biology and Detection ---
    lambda_tilde_i <- as.numeric(w_matrix %*% full_lambda_j)
    psi_i <- 1 - exp(-lambda_tilde_i)
    Z_i <- rbinom(length(psi_i), 1, psi_i)
    site_Z <- setNames(Z_i, rownames(w_matrix))
    
    # Simulate for points directly inside these generated geometries
    final_points$obs_cov1 <- rnorm(nrow(final_points))
    
    # FIX: Filter points to only those that belong to valid sites in the w_matrix
    final_points <- final_points[as.character(final_points$site) %in% rownames(w_matrix), ]
    
    final_points$Z_i <- site_Z[as.character(final_points$site)]
    
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
    
    if (length(valid_M_values) < length(M_values_to_test)) {
      cat(sprintf("  (Note: Available clusters %d limits max M tested)\n", total_sites))
    }
    
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
      
      # 1. Create and specifically align y_wide
      y_df <- occu_df %>% select(site, visit, y) %>%
        pivot_wider(names_from=visit, values_from=y)
      y_df <- y_df[match(rownames(sub_w), y_df$site), ] # CRITICAL: Align rows to w_matrix
      y_wide <- y_df %>% select(-site) %>% as.matrix()
      
      # 2. Create and specifically align obs_cov1_wide
      obs_cov1_df <- occu_df %>% select(site, visit, obs_cov1) %>%
        pivot_wider(names_from=visit, values_from=obs_cov1)
      obs_cov1_df <- obs_cov1_df[match(rownames(sub_w), obs_cov1_df$site), ] # CRITICAL: Align rows
      obs_cov1_wide <- obs_cov1_df %>% select(-site) %>% as.matrix()
      
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
        fm_rep <- try(occuPPM(
          formula = ~obs_cov1 ~ cell_cov1,
          data = umf,
          starts = rand_starts,
          se = TRUE,
          method = selected_optimizer,
          lower = PARAM_LOWER, 
          upper = PARAM_UPPER 
        ), silent = TRUE)
        
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
        
        # We place the cell/pixel level data into our grid
        cell_df <- data.frame(
          x = full_cell_col,
          y = full_cell_row,
          covariate = full_cellCovs$cell_cov1,
          pixel_abundance = full_lambda_j,
          pixel_occupancy = as.factor(full_pixel_Z)
        )
        
        # Geometries for plotting (Unfilled, Red borders only)
        sel_geoms <- final_geoms %>% filter(site %in% sel_sites)
        
        tight_theme <- theme_minimal() + 
          theme(
            axis.title = element_blank(),
            plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt"),
            legend.position = "bottom",
            legend.direction = "horizontal"
          )
        
        p_cov <- ggplot(cell_df, aes(x=x, y=y, fill=covariate)) +
          geom_raster() +
          scale_fill_viridis_c() +
          geom_sf(data=sel_geoms, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          geom_point(data=sf::st_drop_geometry(sub_pts), aes(x=x, y=y), color="black", size=0.5, inherit.aes=FALSE) +
          coord_sf(expand=FALSE, datum=NA) + 
          labs(title=sprintf("Covariate (M=%d)", M_i), fill="Covariate") +
          tight_theme
        
        p_abund <- ggplot(cell_df, aes(x=x, y=y, fill=pixel_abundance)) +
          geom_raster() +
          scale_fill_viridis_c(option = "magma") +
          geom_sf(data=sel_geoms, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          coord_sf(expand=FALSE, datum=NA) + 
          labs(title=sprintf("Abundance (M=%d)", M_i), fill="Abundance") +
          tight_theme
        
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
        theme(legend.position = "bottom", legend.direction = "horizontal") # Add legend to individual plots
    }
    
    p1 <- create_error_plot("beta (state_int)", "State Intercept")
    p2 <- create_error_plot("beta (state_cov1)", "State Slope")
    p3 <- create_error_plot("alpha (det_int)", "Observation Intercept")
    p4 <- create_error_plot("alpha (det_cov1)", "Observation Slope")
    
    # Combine natively without the problematic '&' operator
    combined_error_plot <- (p1 | p2) / (p3 | p4) + 
      patchwork::plot_layout(guides = "collect")
    
    suppressMessages(suppressWarnings({
        ggsave(file.path(output_dir, "error_boxplots.png"), plot = combined_error_plot, dpi = 300, width = 10, height = 10)
    }))
    
    cat("\nSaved results and error boxplots.")
}

cat("\n--- Script Finished Successfully ---\n")