################################################################
# Updated Simulation Experiments: Sampling Extents (eBird Biased)
# Incorporating lat-long, 1to10, and 2to10 methods and GT variants
################################################################

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
library(sqldf)
library(RTriangle)
library(igraph)
library(data.table)
library(ClustGeo)

# Make sure this path correctly points to your dbsc.R helper file
source("R/clustering/dbsc.R")

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Simulation repetitions ---
n_sims <- 10
n_reps <- 30 

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200
# full_grid_dim <- 400

full_n_cells <- full_grid_dim * full_grid_dim # 40000

# --- Fixed Data Point Parameters ---
n_total_points <- 5000
n_single_locs <- 1000
n_double_locs <- 1000
n_hotspot_locs <- 40
hotspot_visits <- 50
hotspot_noise_radius <- 10 # cells

# --- True parameter values ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
true_betas <- c(beta_int = -5.0, beta_cov = 1.0)
# true_betas <- c(beta_int = -3.0, beta_cov = 1.0)

# --- Model settings ---
selected_optimizer <- "nlminb"

PARAM_LOWER <- -10
PARAM_UPPER <- 10
INIT_LOWER <- -2
INIT_UPPER <- 2

# --- Fixed SAC and Skew ---
sac_sigma <- 3 # Medium SAC
n_centers <- 1
centers_scale <- 5
decay_scale <- 30^2


# --- Dynamic Grid Method Settings ---
grid_size <- 5
grid_method_name <- paste0(grid_size, "-cellSq") 

# --- Geometry Buffer Setting ---
buffer_size_cells <- 2

# --- Methods to Run ---
# Moved to the top so you can control which methods run globally
methods_to_run <- c("gt-lat-long", "gt-1to10", "gt-2to10", "lat-long", "1to10", "2to10", grid_method_name, "DBSC", "clustGeo")


# --- 3 Scenarios: Sampling Extents ---
extents <- c("Small" = 1600, "Medium" = 400, "Large" = 100)

# --- ClustGeo Method Settings ---
# Set the kappa percentage (0-100) used for each extent size
kappa_values <- c("Small" = 20, "Medium" = 5, "Large" = 1)

# FORCE absolute path
output_dir <- file.path(getwd(), "output", "simulation_experiments", "updated")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- Simulation Starting ---\n")
cat(sprintf("Running %d simulations for 3 Sampling Extents.\n", n_sims))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", n_sims * length(extents) * 6 * n_reps)) # 6 methods now

##########
# 3. Helpers for Simulation Geometry (2 cell buffer)
##########

sim_create_geometries <- function(obs_df, r_trend, buffer_cells = 2) {
  pts <- sf::st_as_sf(obs_df, coords = c("reported_x", "reported_y"))
  bbox_poly <- sf::st_as_sfc(sf::st_bbox(pts) + buffer_cells * 3)
  
  v_res <- try(sf::st_voronoi(sf::st_union(pts), envelope = bbox_poly), silent=TRUE)
  if(inherits(v_res, "try-error")) return(data.frame())
  
  voronoi_tiles <- sf::st_collection_extract(v_res, "POLYGON") %>% sf::st_sf()
  
  voronoi_w_id <- sf::st_join(voronoi_tiles, pts, join = sf::st_contains)
  site_territories <- voronoi_w_id %>% dplyr::group_by(site) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>% sf::st_make_valid()
  
  site_buffers <- pts %>% dplyr::group_by(site) %>%
    dplyr::summarise(geometry = sf::st_convex_hull(sf::st_union(geometry)), .groups = "drop") %>%
    sf::st_buffer(dist = buffer_cells) %>% sf::st_make_valid()
  
  site_territories <- dplyr::rename(site_territories, site_t = site)
  site_buffers <- dplyr::rename(site_buffers, site_b = site)
  
  sf::st_agr(site_territories) <- "constant"
  sf::st_agr(site_buffers) <- "constant"
  
  intersections <- suppressWarnings(sf::st_intersection(site_territories, site_buffers))
  final_geoms <- intersections %>%
    dplyr::filter(site_t == site_b) %>%
    dplyr::rename(site = site_t) %>%
    dplyr::select(site, geometry) %>%
    sf::st_make_valid() %>%
    sf::st_simplify(dTolerance = 0.1, preserveTopology = TRUE)
  
  geom_type <- sf::st_geometry_type(final_geoms, by_geometry = FALSE)
  if (inherits(geom_type, "GEOMETRYCOLLECTION") || any(grepl("COLLECTION", geom_type))) {
    final_geoms <- sf::st_collection_extract(final_geoms, "POLYGON")
  }
  return(final_geoms)
}

sim_disjoint <- function(site_geoms_sf, point_data_df) {
  sites_split <- site_geoms_sf %>%
    sf::st_make_valid() %>%
    suppressWarnings(sf::st_cast("MULTIPOLYGON")) %>% 
    suppressWarnings(sf::st_cast("POLYGON"))
  
  if (nrow(sites_split) == nrow(site_geoms_sf)) {
    return(list(geoms = site_geoms_sf, data = point_data_df))
  }
  
  sites_split <- sites_split %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(sub_id = dplyr::row_number(), new_site_id = paste0(site, "_", sub_id)) %>%
    dplyr::ungroup()
  
  points_sf <- sf::st_as_sf(point_data_df, coords = c("reported_x", "reported_y"))
  join_res <- suppressWarnings(sf::st_join(points_sf, sites_split["new_site_id"], join = sf::st_intersects, left = TRUE))
  join_res <- join_res[!duplicated(join_res$checklist_id), ]
  
  point_data_updated <- point_data_df
  match_idx <- match(point_data_updated$checklist_id, join_res$checklist_id)
  new_ids <- join_res$new_site_id[match_idx]
  
  point_data_updated$site <- ifelse(!is.na(new_ids), new_ids, as.character(point_data_updated$site))
  
  occupied_new_ids <- unique(na.omit(new_ids))
  sites_final <- sites_split %>%
    dplyr::filter(new_site_id %in% occupied_new_ids) %>%
    dplyr::select(-site, -sub_id) %>%
    dplyr::rename(site = new_site_id) %>%
    dplyr::select(site, geometry)
  
  return(list(geoms = sites_final, data = point_data_updated))
}

sim_overlap <- function(site_geoms_sf, r_trend) {
  site_vect <- terra::vect(site_geoms_sf)
  overlap_df <- terra::extract(r_trend, site_vect, cells = TRUE, exact = TRUE, ID = TRUE)
  overlap_df <- overlap_df[!is.na(overlap_df[[names(r_trend)[1]]]), ]
  
  if (!("fraction" %in% names(overlap_df))) overlap_df$fraction <- 1.0
  overlap_df$w_area <- overlap_df$fraction * 1
  
  n_sites <- nrow(site_geoms_sf)
  n_cells <- terra::ncell(r_trend)
  
  w <- Matrix::sparseMatrix(
    i = overlap_df$ID,
    j = overlap_df$cell,
    x = overlap_df$w_area,
    dims = c(n_sites, n_cells),
    dimnames = list(site_geoms_sf$site, NULL)
  )
  return(w)
}

fit_clustered_model <- function(clustered_obs, w_matrix, full_cellCovs) {
  site_counts <- table(clustered_obs$site)
  active_sites <- as.character(names(site_counts))
  M_active <- length(active_sites)
  if (M_active == 0) return(NULL)
  J_max <- max(site_counts)
  
  y_mat <- matrix(NA, M_active, J_max)
  obs_mat <- matrix(NA, M_active, J_max)
  
  for(i in 1:M_active) {
    s_idx <- active_sites[i]
    s_obs <- clustered_obs[clustered_obs$site == s_idx, ]
    J_i <- nrow(s_obs)
    y_mat[i, 1:J_i] <- s_obs$detection
    obs_mat[i, 1:J_i] <- s_obs$obs_cov
  }
  
  w_active <- w_matrix[active_sites, , drop = FALSE]
  
  umf <- unmarkedFrameOccuPPM(
    y = y_mat, obsCovs = list(obs_cov1 = obs_mat), 
    cellCovs = full_cellCovs, w = w_active 
  )
  
  best_fm <- NULL
  min_nll <- Inf
  n_params <- length(true_alphas) + length(true_betas)
  
  for (rep in 1:n_reps) {
    rand_starts <- runif(n_params, min = INIT_LOWER, max = INIT_UPPER)
    fm_rep <- try(occuPPM(
      formula = ~obs_cov1 ~ cell_cov1,
      data = umf, starts = rand_starts, se = TRUE,
      method = selected_optimizer, lower = PARAM_LOWER, upper = PARAM_UPPER
    ), silent = TRUE)
    
    if (!inherits(fm_rep, "try-error")) {
      if (fm_rep@negLogLike < min_nll) {
        min_nll <- fm_rep@negLogLike
        best_fm <- fm_rep
      }
    }
  }
  return(best_fm)
}


##########
# 4. Pre-generate Coordinates & Seeds
##########

full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1
cell_coords <- data.frame(x = full_cell_col, y = full_cell_row)

set.seed(123)
cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(n_centers, 0, full_grid_dim), 
    y = runif(n_centers, 0, full_grid_dim)
  )
}

##########
# 5. Main Simulation Loop
##########

all_results <- list()
plot_data <- list()
geom_plot_data <- list()

for (sim in 1:n_sims) {
  cat(sprintf("\n=== Sim %d of %d ===\n", sim, n_sims))
  
  r_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                         xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim,
                         vals=rnorm(full_n_cells))
  terra::crs(r_noise) <- ""
  
  fw <- terra::focalMat(r_noise, sac_sigma, type = "Gauss")
  r_smooth <- terra::focal(r_noise, w = fw, fun = sum, na.rm = TRUE)
  terra::values(r_smooth) <- as.vector(scale(terra::values(r_smooth)))
  
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
  
  full_cellCovs <- data.frame(cell_cov1 = terra::values(r_trend, mat=FALSE))
  if(any(is.na(full_cellCovs$cell_cov1))) full_cellCovs$cell_cov1[is.na(full_cellCovs$cell_cov1)] <- 0
  
  full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
  full_lambda_j <- exp(full_X_cell %*% true_betas)
  
  full_N_j <- rpois(full_n_cells, full_lambda_j)
  full_Z_j <- factor(ifelse(full_N_j > 0, 1, 0), levels = c("0", "1"))
  
  r_acc_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                             xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim,
                             vals=rnorm(full_n_cells))
  r_acc <- terra::focal(r_acc_noise, w = fw, fun = sum, na.rm = TRUE)
  acc_vals <- as.vector(scale(terra::values(r_acc)))
  acc_weights <- exp(acc_vals * 2) 
  
  set.seed(100 + sim)
  sampled_loc_ids <- sample(1:full_n_cells, n_hotspot_locs + n_double_locs + n_single_locs, prob = acc_weights)
  
  hotspot_idx <- sampled_loc_ids[1:n_hotspot_locs]
  double_idx <- sampled_loc_ids[(n_hotspot_locs+1):(n_hotspot_locs+n_double_locs)]
  single_idx <- sampled_loc_ids[(n_hotspot_locs+n_double_locs+1):length(sampled_loc_ids)]
  
  obs_list <- list()
  obs_list[[1]] <- data.frame(reported_cell = single_idx, true_cell = single_idx)
  obs_list[[2]] <- data.frame(reported_cell = rep(double_idx, each=2), true_cell = rep(double_idx, each=2))
  
  hotspot_reported <- rep(hotspot_idx, each = hotspot_visits)
  hotspot_true <- numeric(length(hotspot_reported))
  
  for(i in seq_along(hotspot_idx)) {
    base_idx <- hotspot_idx[i]
    base_x <- full_cell_col[base_idx]
    base_y <- full_cell_row[base_idx]
    
    t_cells <- numeric(hotspot_visits)
    t_cells[1:3] <- base_idx
    for(j in 4:hotspot_visits) {
      angle <- runif(1, 0, 2*pi)
      r <- sqrt(runif(1, 0, hotspot_noise_radius^2))
      nx <- round(base_x + r * cos(angle))
      ny <- round(base_y + r * sin(angle))
      nx <- max(1, min(full_grid_dim, nx))
      ny <- max(1, min(full_grid_dim, ny))
      t_cells[j] <- (ny - 1) * full_grid_dim + nx
    }
    hotspot_true[((i-1)*hotspot_visits + 1):(i*hotspot_visits)] <- t_cells
  }
  obs_list[[3]] <- data.frame(reported_cell = hotspot_reported, true_cell = hotspot_true)
  
  all_obs <- do.call(rbind, obs_list)
  all_obs$reported_x <- full_cell_col[all_obs$reported_cell]
  all_obs$reported_y <- full_cell_row[all_obs$reported_cell]
  all_obs$checklist_id <- 1:nrow(all_obs)
  obs_sf <- sf::st_as_sf(all_obs, coords = c("reported_x", "reported_y"))
  
  if (sim == 1) {
    plot_data[["Cell"]] <- data.frame(
      x = full_cell_col, y = full_cell_row, covariate = full_cellCovs$cell_cov1,
      abundance = full_N_j, occupancy = full_Z_j
    )
  }
  
  if (sim == 1) cat("Generating Density-Weighted Voronoi Site Geometries...\n")
  site_definitions <- list()
  prob_weights <- exp(full_cellCovs$cell_cov1 * 0.5)
  
  for (ext_name in names(extents)) {
    K <- extents[[ext_name]]
    seed_idx <- sample(1:full_n_cells, K, prob = prob_weights)
    seed_pts <- cell_coords[seed_idx, ]
    
    site_ids <- rep(1, full_n_cells)
    min_dists <- rep(Inf, full_n_cells)
    for (k in 1:K) {
      dists_sq <- (cell_coords$x - seed_pts$x[k])^2 + (cell_coords$y - seed_pts$y[k])^2
      update_idx <- dists_sq < min_dists
      site_ids[update_idx] <- k
      min_dists[update_idx] <- dists_sq[update_idx]
    }
    
    w <- Matrix::sparseMatrix(
      i = site_ids, j = 1:full_n_cells, x = 1, dims = c(K, full_n_cells)
    )
    
    site_sf <- NULL
    if (sim == 1) {
      xyz_df <- data.frame(x = full_cell_col, y = full_cell_row, z = site_ids)
      r_sites <- terra::rast(xyz_df, type = "xyz")
      site_polys <- terra::as.polygons(r_sites, dissolve = TRUE)
      site_sf <- sf::st_as_sf(site_polys)
      colnames(site_sf)[1] <- "site"
    }
    
    site_definitions[[ext_name]] <- list(K = K, site_ids = site_ids, w = w, site_sf = site_sf)
  }
  
  for (ext_name in names(extents)) {
    cat(sprintf("  - Extent: %s\n", ext_name))
    
    def <- site_definitions[[ext_name]]
    M <- def$K
    w <- def$w
    rownames(w) <- as.character(1:M)
    
    lambda_tilde_i <- as.numeric(w %*% full_lambda_j)
    N_i <- rpois(M, lambda_tilde_i)
    Z_i <- factor(ifelse(N_i > 0, 1, 0), levels = c("0", "1"))
    
    all_obs$reported_site <- def$site_ids[all_obs$reported_cell]
    all_obs$true_site <- def$site_ids[all_obs$true_cell]
    
    # 1. GENERATE FIXED DETECTIONS BASED ON GROUND TRUTH POLYGONS
    all_obs$obs_cov <- rnorm(nrow(all_obs))
    all_obs$detection <- 0
    for(r_idx in 1:nrow(all_obs)) {
      Z_t <- Z_i[all_obs$true_site[r_idx]]
      if (Z_t == "1") {
        logit_p <- true_alphas[1] + true_alphas[2] * all_obs$obs_cov[r_idx]
        all_obs$detection[r_idx] <- rbinom(1, 1, plogis(logit_p))
      }
    }
    
    if (sim == 1) {
      active_sites <- unique(all_obs$reported_site)
      agg_cov <- as.numeric((w %*% full_cellCovs$cell_cov1) / rowSums(w))
      
      sf_data <- def$site_sf
      sf_data <- sf_data[order(sf_data$site), ] 
      
      sf_data$is_active <- (1:M %in% active_sites)
      sf_data$covariate <- agg_cov
      sf_data$covariate[!sf_data$is_active] <- NA
      sf_data$abundance <- N_i
      sf_data$abundance[!sf_data$is_active] <- NA
      sf_data$occupancy <- Z_i
      sf_data$occupancy[!sf_data$is_active] <- NA
      
      plot_data[[ext_name]] <- list(sf_data = sf_data, obs_sf = obs_sf)
    }

    # 2. RUN ALL METHODS FOR THIS EXTENT
    
    for (m_name in methods_to_run) {
      cat(sprintf("    * Method: %s\n", m_name))
      
      # Step A: Filter observations based on method logic
      if (grepl("lat-long", m_name)) {
        obs_m <- all_obs
      } else if (grepl("1to10", m_name)) {
        obs_m <- all_obs %>% dplyr::group_by(reported_x, reported_y) %>% dplyr::slice_sample(n = 10) %>% dplyr::ungroup() %>% as.data.frame()
      } else if (grepl("2to10", m_name) || m_name == grid_method_name || m_name == "DBSC" || m_name == "clustGeo") {
        # 2to10, dynamic grid, DBSC, and clustGeo all start with the 2-to-10 filtering baseline
        obs_m <- all_obs %>% dplyr::group_by(reported_x, reported_y) %>% dplyr::filter(n() >= 2) %>% dplyr::slice_sample(n = 10) %>% dplyr::ungroup() %>% as.data.frame()
      }
      
      if (nrow(obs_m) == 0) next
      
      # Step B: Assign geometry and compute parameters
      if (startsWith(m_name, "gt-")) {
        # Ground-truth methods
        obs_m$site <- as.character(obs_m$reported_site)
        fm_m <- fit_clustered_model(obs_m, w, full_cellCovs)
        
      } else if (m_name == grid_method_name) {
        # Dynamic Cell Grid Method 
        obs_m$grid_x <- (obs_m$reported_x - 1) %/% grid_size
        obs_m$grid_y <- (obs_m$reported_y - 1) %/% grid_size
        obs_m$site <- paste0("grid_", obs_m$grid_x, "_", obs_m$grid_y)
        
        # Build exact square geometries for occupied sites
        active_sites <- unique(obs_m$site)
        grid_polys <- list()
        for (s in active_sites) {
          parts <- strsplit(gsub("grid_", "", s), "_")[[1]]
          gx <- as.numeric(parts[1])
          gy <- as.numeric(parts[2])
          # Grid bounds aligned with cell centers
          xmin <- gx * grid_size + 0.5; xmax <- xmin + grid_size
          ymin <- gy * grid_size + 0.5; ymax <- ymin + grid_size
          poly <- sf::st_polygon(list(matrix(c(xmin,ymin, xmax,ymin, xmax,ymax, xmin,ymax, xmin,ymin), ncol=2, byrow=TRUE)))
          grid_polys[[s]] <- sf::st_sf(site = s, geometry = sf::st_sfc(poly))
        }
        geoms <- do.call(rbind, grid_polys)
        
        # Standard processing 
        split_res <- sim_disjoint(geoms, obs_m)
        w_m <- sim_overlap(split_res$geoms, r_trend)
        fm_m <- fit_clustered_model(split_res$data, w_m, full_cellCovs)
        
      } else if (m_name == "DBSC") {
        # DBSC Method 
        # Spoof names so runDBSC() works seamlessly
        obs_m$latitude <- obs_m$reported_y
        obs_m$longitude <- obs_m$reported_x
        obs_m$cell_cov1 <- full_cellCovs$cell_cov1[obs_m$reported_cell]
        
        # Run clustering
        obs_m <- runDBSC(obs_m, "cell_cov1")
        obs_m$site <- as.character(obs_m$site)
        
        # Apply standard Polygonization (Voronoi + Convex Hull + Cell Buffer)
        geoms <- sim_create_geometries(obs_m, r_trend, buffer_cells = buffer_size_cells)
        if (nrow(geoms) == 0) next
        
        split_res <- sim_disjoint(geoms, obs_m)
        w_m <- sim_overlap(split_res$geoms, r_trend)
        fm_m <- fit_clustered_model(split_res$data, w_m, full_cellCovs)
        
      } else if (m_name == "clustGeo") {
        # clustGeo Method 
        
        # Set kappa dynamically based on the current extent
        kappa <- kappa_values[[ext_name]]
        # Prepare distinct locations and covariates
        obs_m$cell_cov1 <- full_cellCovs$cell_cov1[obs_m$reported_cell]
        obs_locs <- obs_m %>% dplyr::distinct(reported_x, reported_y, .keep_all = TRUE)
        
        # Create distance matrices
        env_dist <- dist(scale(obs_locs$cell_cov1))
        
        df_for_dist <- obs_locs
        df_for_dist$y_scaled <- as.numeric(scale(df_for_dist$reported_y))
        df_for_dist$x_scaled <- as.numeric(scale(df_for_dist$reported_x))
        geo_dist <- dist(df_for_dist[, c("y_scaled", "x_scaled")])
        
        # Run clustering
        tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = 0.5)
        percent <- kappa / 100.0
        K <- max(2, round(nrow(obs_locs) * percent))
        obs_locs$site <- cutree(tree, K)
        
        # Map sites back to observations
        obs_m <- obs_m %>% dplyr::left_join(obs_locs[, c("reported_x", "reported_y", "site")], by = c("reported_x", "reported_y"))
        obs_m$site <- as.character(obs_m$site)
        
        # Apply standard Polygonization (Voronoi + Convex Hull + Cell Buffer)
        geoms <- sim_create_geometries(obs_m, r_trend, buffer_cells = buffer_size_cells)
        if (nrow(geoms) == 0) next
        
        split_res <- sim_disjoint(geoms, obs_m)
        w_m <- sim_overlap(split_res$geoms, r_trend)
        fm_m <- fit_clustered_model(split_res$data, w_m, full_cellCovs)
      }
      else {
        # Standard clustering methods (lat-long, 1to10, 2to10)
        obs_m$site <- paste0(obs_m$reported_x, "_", obs_m$reported_y)
        
        geoms <- sim_create_geometries(obs_m, r_trend, buffer_cells = buffer_size_cells)
        if (nrow(geoms) == 0) next
        
        split_res <- sim_disjoint(geoms, obs_m)
        w_m <- sim_overlap(split_res$geoms, r_trend)
        fm_m <- fit_clustered_model(split_res$data, w_m, full_cellCovs)
      }
      
      if (sim == 1) {
        plot_key <- NULL
        # Capture fixed methods once during the "Small" extent loop to avoid overwriting
        if (m_name == grid_method_name && ext_name == "Small") plot_key <- grid_method_name
        if (m_name == "DBSC" && ext_name == "Small") plot_key <- "DBSC"
        if (m_name == "lat-long" && ext_name == "Small") plot_key <- "lat-long"
        if (m_name == "1to10" && ext_name == "Small") plot_key <- "1to10"
        if (m_name == "2to10" && ext_name == "Small") plot_key <- "2to10"
        
        # Capture clustGeo dynamically across extents to get all 3 kappas
        if (m_name == "clustGeo") {
          plot_key <- paste0("clustGeo-0.5-", kappa_values[[ext_name]])
        }

        if (!is.null(plot_key) && exists("split_res")) {
           geom_plot_data[[plot_key]] <- list(
               geoms = split_res$geoms,
               pts = sf::st_as_sf(obs_m, coords = c("reported_x", "reported_y"))
           )
        }
      }

      est_val <- if(is.null(fm_m)) c(NA,NA,NA,NA) else c(coef(fm_m, 'det'), coef(fm_m, 'state'))
      
      all_results[[length(all_results) + 1]] <- data.frame(
        Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
        True_Value = c(true_alphas, true_betas),
        Estimated_Value = est_val,
        Extent = ext_name,
        Method = m_name,
        sim_id = sim
      )
    }
  }
}

##########
# 6. Generate and Save 4x3 Plot (Sampling Extents) - Ground-Truth Visualization
##########
cat("\nGenerating 4x3 Spatial Plot (sampling_extents.png)...\n")

cov_limits <- range(c(plot_data$Cell$covariate, 
                      plot_data$Small$sf_data$covariate, 
                      plot_data$Medium$sf_data$covariate, 
                      plot_data$Large$sf_data$covariate), na.rm=TRUE)

abund_limits <- range(c(plot_data$Cell$abundance, 
                        plot_data$Small$sf_data$abundance, 
                        plot_data$Medium$sf_data$abundance, 
                        plot_data$Large$sf_data$abundance), na.rm=TRUE)

ggplot2::theme_set(ggplot2::theme_minimal() + ggplot2::theme(
  legend.position = "bottom",
  legend.justification = "center",
  legend.box = "horizontal",
  legend.box.just = "center",
  legend.spacing.x = ggplot2::unit(2.5, "cm"), 
  legend.margin = ggplot2::margin(t = 5, l = 15),
  legend.title.align = 0.5, 
  legend.title = ggplot2::element_text(size=14),
  legend.text = ggplot2::element_text(size=12)
))

base_theme <- ggplot2::theme(
  axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
  panel.grid = ggplot2::element_blank(), plot.margin = ggplot2::margin(2, 2, 2, 2, "pt"),
  plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold"),
  axis.title.y = ggplot2::element_text(size = 14, face = "bold", angle = 90, vjust = 0.5)
)

guide_cont <- ggplot2::guide_colorbar(
  direction = "horizontal", title.position = "top", title.hjust = 0.5,
  barwidth = ggplot2::unit(5.0, "cm"), barheight = ggplot2::unit(0.6, "cm")
)

guide_disc <- ggplot2::guide_legend(
  direction = "horizontal", title.position = "top", title.hjust = 0.5,
  keywidth = ggplot2::unit(1.0, "cm"), keyheight = ggplot2::unit(0.6, "cm")
)

build_row <- function(data_name, row_title, show_titles=FALSE) {
  if (data_name == "Cell") {
    df <- plot_data[[data_name]]
    p1 <- ggplot(df, aes(x=x, y=y, fill=covariate)) + geom_raster() + scale_fill_viridis_c(name="Covariate", limits=cov_limits, guide=guide_cont) + base_theme + labs(y = row_title, x = NULL) + coord_fixed(expand=FALSE)
    p2 <- ggplot(df, aes(x=x, y=y, fill=abundance)) + geom_raster() + scale_fill_viridis_c(option="magma", name="Abundance", limits=abund_limits, guide=guide_cont) + base_theme + labs(y = NULL, x = NULL) + coord_fixed(expand=FALSE)
    p3 <- ggplot(df, aes(x=x, y=y, fill=occupancy)) + geom_raster() + scale_fill_manual(values=c("0"="#440154FF", "1"="#FDE725FF"), name="           Occupancy           ", drop=FALSE, guide=guide_disc) + base_theme + labs(y = NULL, x = NULL) + coord_fixed(expand=FALSE)
  } else {
    df <- plot_data[[data_name]]$sf_data
    obs_pts <- plot_data[[data_name]]$obs_sf
    p1 <- ggplot(df) + geom_sf(aes(fill=covariate), color="black", linewidth=0.1) + geom_sf(data=obs_pts, color="red", size=0.2, alpha=0.5) + scale_fill_viridis_c(name="Covariate", limits=cov_limits, guide=guide_cont, na.value="white") + base_theme + labs(y = row_title, x = NULL) + coord_sf(expand=FALSE)
    p2 <- ggplot(df) + geom_sf(aes(fill=abundance), color="black", linewidth=0.1) + scale_fill_viridis_c(option="magma", name="Abundance", limits=abund_limits, guide=guide_cont, na.value="white") + base_theme + labs(y = NULL, x = NULL) + coord_sf(expand=FALSE)
    p3 <- ggplot(df) + geom_sf(aes(fill=occupancy), color="black", linewidth=0.1, show.legend=FALSE) + scale_fill_manual(values=c("0"="#440154FF", "1"="#FDE725FF"), name="           Occupancy           ", drop=FALSE, guide=guide_disc, na.value="white") + base_theme + labs(y = NULL, x = NULL) + coord_sf(expand=FALSE)
  }
  
  if(show_titles) {
    p1 <- p1 + ggtitle("Covariate")
    p2 <- p2 + ggtitle("Abundance")
    p3 <- p3 + ggtitle("Occupancy")
  }
  return(list(p1, p2, p3))
}

row1 <- build_row("Cell", "Simulated Species", show_titles=TRUE)
row2 <- build_row("Small", "Small Sampling Extents")
row3 <- build_row("Medium", "Medium Sampling Extents")
row4 <- build_row("Large", "Large Sampling Extents")

comb_plot <- patchwork::wrap_plots(c(row1, row2, row3, row4), ncol=3) + patchwork::plot_layout(guides="collect")
ggsave(file.path(output_dir, "sampling_extents.png"), plot=comb_plot, width=10, height=14, dpi=300)

##########
# 6.5 Generate 4x3 Plot (Site Geometries)
##########
cat("Generating 4x3 Site Geometries Plot (site_geometries.png)...\n")

# Helper to ensure standard plotting aesthetics across all 12 panels
plot_geom_panel <- function(geoms_sf, pts_sf, title) {
  p <- ggplot()
  if (!is.null(geoms_sf)) {
    p <- p + geom_sf(data = geoms_sf, fill = "lightblue", alpha = 0.3, color = "black", linewidth = 0.2)
  }
  if (!is.null(pts_sf)) {
    p <- p + geom_sf(data = pts_sf, color = "red", size = 0.4, alpha = 0.6)
  }

  p <- p +
    base_theme +
    ggtitle(title) +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(t = 10, b = 5)),
          plot.margin = margin(t = 10, r = 5, b = 5, l = 5), # Increases space above the plot
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    coord_sf(xlim = c(0, full_grid_dim), ylim = c(0, full_grid_dim), expand = FALSE)
  return(p)
}

# Row 1: All Points | Grid Method | DBSC
p_r1c1 <- plot_geom_panel(NULL, plot_data$Small$obs_sf, "Observation Points")
p_r1c2 <- plot_geom_panel(geom_plot_data[[grid_method_name]]$geoms, geom_plot_data[[grid_method_name]]$pts, paste0(grid_size, "-cellSq"))
p_r1c3 <- plot_geom_panel(geom_plot_data[["DBSC"]]$geoms, geom_plot_data[["DBSC"]]$pts, "DBSC")

# Row 2: Standard Filtering (lat-long, 1to10, 2to10)
p_r2c1 <- plot_geom_panel(geom_plot_data[["lat-long"]]$geoms, geom_plot_data[["lat-long"]]$pts, "lat-long")
p_r2c2 <- plot_geom_panel(geom_plot_data[["1to10"]]$geoms, geom_plot_data[["1to10"]]$pts, "1to10")
p_r2c3 <- plot_geom_panel(geom_plot_data[["2to10"]]$geoms, geom_plot_data[["2to10"]]$pts, "2to10")

# Row 3: ClustGeo at different extents/kappas (Dynamically pulled from kappa_values)
cg_key_small <- paste0("clustGeo-0.5-", kappa_values[["Small"]])
cg_key_med   <- paste0("clustGeo-0.5-", kappa_values[["Medium"]])
cg_key_large <- paste0("clustGeo-0.5-", kappa_values[["Large"]])

p_r3c1 <- plot_geom_panel(geom_plot_data[[cg_key_small]]$geoms, geom_plot_data[[cg_key_small]]$pts, cg_key_small)
p_r3c2 <- plot_geom_panel(geom_plot_data[[cg_key_med]]$geoms, geom_plot_data[[cg_key_med]]$pts, cg_key_med)
p_r3c3 <- plot_geom_panel(geom_plot_data[[cg_key_large]]$geoms, geom_plot_data[[cg_key_large]]$pts, cg_key_large)

# Row 4: Sampling Extents (Using ALL points)
p_r4c1 <- plot_geom_panel(plot_data$Small$sf_data[plot_data$Small$sf_data$is_active, ], plot_data$Small$obs_sf, "Small Sampling Extents")
p_r4c2 <- plot_geom_panel(plot_data$Medium$sf_data[plot_data$Medium$sf_data$is_active, ], plot_data$Medium$obs_sf, "Medium Sampling Extents")
p_r4c3 <- plot_geom_panel(plot_data$Large$sf_data[plot_data$Large$sf_data$is_active, ], plot_data$Large$obs_sf, "Large Sampling Extents")

# Combine and save
comb_geoms <- (p_r1c1 | p_r1c2 | p_r1c3) /
              (p_r2c1 | p_r2c2 | p_r2c3) /
              (p_r3c1 | p_r3c2 | p_r3c3) /
              (p_r4c1 | p_r4c2 | p_r4c3)

ggsave(file.path(output_dir, "site_geometries.png"), plot = comb_geoms, width = 10, height = 16, dpi = 300)


##########
# 7. Process Results & Generate Error Boxplots
##########
cat("Saving Results & Generating Error Boxplots...\n")

res_df <- do.call(rbind, all_results)
res_df <- res_df[!is.na(res_df$Estimated_Value), ]

res_df$Error <- res_df$True_Value - res_df$Estimated_Value
res_df$Extent <- factor(res_df$Extent, levels = c("Small", "Medium", "Large"))

# Assign proper levels to match generation
res_df$Method <- factor(res_df$Method, levels = c("gt-lat-long", "gt-1to10", "gt-2to10", "lat-long", "1to10", "2to10", grid_method_name, "DBSC", "clustGeo"))

write.csv(res_df, file.path(output_dir, "params_updated.csv"), row.names = FALSE)


# ----------------------------------------------------------------------
# PLOT 1: error_boxplots.png (Ground-truth variants isolated)
# ----------------------------------------------------------------------
res_df_gt <- res_df[startsWith(as.character(res_df$Method), "gt-"), ]

# Remove the "gt-" prefix to create clean labels for the legend
res_df_gt$Method_Label <- factor(
  gsub("gt-", "", as.character(res_df_gt$Method)),
  levels = c("lat-long", "1to10", "2to10")
)

gt_method_colors <- c("lat-long" = "navy", "1to10" = "cyan", "2to10" = "pink")

create_error_plot_gt <- function(param_name, title) {
  ggplot(res_df_gt[res_df_gt$Parameter == param_name, ], 
         # Swapped x and fill here:
         aes(x = Extent, y = Error, fill = Method_Label)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    # Updated to the new colors and added a legend title:
    scale_fill_manual(name = "Data Points", values = gt_method_colors) +
    # Updated the x-axis label:
    labs(title = title, x = "Sampling Extent", y = "Error (True - Estimate)") +
    theme_bw() + theme(legend.position = "bottom")
}

p_beta0_gt <- create_error_plot_gt("beta (state_int)", "State Intercept")
p_beta1_gt <- create_error_plot_gt("beta (state_cov1)", "State Slope")
p_alpha0_gt <- create_error_plot_gt("alpha (det_int)", "Observation Intercept")
p_alpha1_gt <- create_error_plot_gt("alpha (det_cov1)", "Observation Slope")

combined_error_plot_gt <- (p_beta0_gt | p_beta1_gt) / (p_alpha0_gt | p_alpha1_gt) + 
  patchwork::plot_layout(guides = "collect")

ggsave(file.path(output_dir, "error_boxplots.png"), plot = combined_error_plot_gt, dpi = 300, width = 10, height = 10)

# ----------------------------------------------------------------------
# PLOT 2: error_boxplots_all.png (reference-2to10 vs clustered)
# ----------------------------------------------------------------------
# Filter only to the methods required for this plot
res_df_all <- res_df[res_df$Method %in% c("gt-2to10", "lat-long", "1to10", "2to10", grid_method_name, "DBSC", "clustGeo"), ]

# Map gt-2to10 to the requested display name
res_df_all$Method_Label <- as.character(res_df_all$Method)
res_df_all$Method_Label[res_df_all$Method_Label == "gt-2to10"] <- "reference-2to10"

# Enforce factor levels for ordering in the plot/legend
res_df_all$Method_Label <- factor(
  res_df_all$Method_Label, 
  levels = c("reference-2to10", "lat-long", "1to10", "2to10", grid_method_name, "DBSC", "clustGeo")
)

# Apply requested color scheme 
method_colors <- c(
  "reference-2to10" = "red", 
  "lat-long" = "navy", 
  "1to10" = "cyan", 
  "2to10" = "pink",
  "DBSC" = "darkgrey",
  "clustGeo" = "forestgreen"
)
# Add the dynamic grid color to the named vector
method_colors[[grid_method_name]] <- "orange"

create_error_plot_all <- function(param_name, title) {
  ggplot(res_df_all[res_df_all$Parameter == param_name, ], 
         aes(x = Extent, y = Error, fill = Method_Label)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(name = "Method", values = method_colors) +
    labs(title = title, x = "Sampling Extent", y = "Error (True - Estimate)") +
    theme_bw() + theme(legend.position = "bottom")
}

p_beta0_all <- create_error_plot_all("beta (state_int)", "State Intercept")
p_beta1_all <- create_error_plot_all("beta (state_cov1)", "State Slope")
p_alpha0_all <- create_error_plot_all("alpha (det_int)", "Observation Intercept")
p_alpha1_all <- create_error_plot_all("alpha (det_cov1)", "Observation Slope")

combined_error_plot_all <- (p_beta0_all | p_beta1_all) / (p_alpha0_all | p_alpha1_all) + 
  patchwork::plot_layout(guides = "collect")

ggsave(file.path(output_dir, "error_boxplots_all.png"), plot = combined_error_plot_all, dpi = 300, width = 12, height = 10)

cat("--- Script Finished Successfully ---\n")