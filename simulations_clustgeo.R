################################################################
# Simulation Experiments using clustGeo Site Definition
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
library(terra) 
library(Matrix)
library(patchwork)
library(scales) 
library(sf)
library(ClustGeo)

# Source required helpers
source(file.path("R", "utils.R"))
source(file.path("R", "clustering_helpers.R"))
source(file.path("R", "model_helpers.R"))

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Simulation repetitions ---
n_sims <- 3

# --- Model fitting repetitions ---
n_reps <- 30 

# --- Point Generation & ClustGeo Parameters ---
max_N_value <- 10000
kappa_for_clustgeo <- 10    # 10% of 10000 = 1000 initial sites
alpha_clustgeo <- 0.5

# --- Spatial & Buffer Parameters ---
full_grid_dim <- 200 
cell_size <- 100 
max_coord <- full_grid_dim * cell_size # <-- MISSING VARIABLE ADDED HERE
buffer_pixels <- 2
buffer_m <- buffer_pixels * cell_size

# Base Raster CRS (Albers to satisfy both terra and sf geometry processing)
albers_crs <- "+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

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
n_centers <- 1
centers_scale <- 5
decay_scale <- (30 * cell_size)^2

output_dir <- file.path("output", "simulation_experiments", "clustgeo")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- ClustGeo Simulation Starting ---\n")
cat(sprintf("Points generated per sim: %d\n", max_N_value))
cat(sprintf("ClustGeo Kappa: %d\n", kappa_for_clustgeo))
cat(sprintf("Buffer size: %d pixels (%dm)\n", buffer_pixels, buffer_m))
cat(sprintf("Running %d full simulations per SAC level.\n", n_sims))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", length(sac_levels) * n_sims * length(M_values_to_test) * n_reps))


##########
# 3. Pre-generate Fixed Seeds
##########

cat("Pre-generating fixed COVARIATE skew seeds...\n")
cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(n_centers, 0, max_coord), 
    y = runif(n_centers, 0, max_coord)
  )
}

##########
# 4. Initialize Storage
##########

results_storage <- list()

##########
# 5. OUTER LOOPS: SAC -> SIM
##########

for (sac_level in sac_levels) {
  current_sigma <- sac_sigmas[[sac_level]]
  
  for (sim in 1:n_sims) {
    cat(sprintf("\n=== SAC: %s | Sim %d of %d ===\n", sac_level, sim, n_sims))
    
    ##########
    # 6. Generate Landscape Data (Cells)
    ##########
    cat("  Generating landscape & raster covariates...\n")
    
    # Generate Base Spatial Noise
    r_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                           xmin=0, xmax=max_coord, ymin=0, ymax=max_coord,
                           vals=rnorm(full_grid_dim^2), crs=albers_crs)
    names(r_noise) <- "cell_cov1"
    
    if (current_sigma > 0) {
      # Multiply sigma by cell_size so 3 becomes 300m (i.e., 3 pixels)
      fw <- terra::focalMat(r_noise, current_sigma * cell_size, type = "Gauss")
      r_smooth <- terra::focal(r_noise, w = fw, fun = sum, na.rm = TRUE)
    } else {
      r_smooth <- r_noise
    }
    terra::values(r_smooth) <- as.vector(scale(terra::values(r_smooth)))
    
    # Apply Skew / Trend
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
    
    # Calculate global lambda_j
    full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
    full_lambda_j <- exp(full_X_cell %*% true_betas)
    
    ##########
    # 7. Generate Points & Extract Covariates
    ##########
    cat("  Generating continuous points and applying clustGeo...\n")
    
    points_df <- data.frame(
      checklist_id = 1:max_N_value,
      x = runif(max_N_value, 0, max_coord),
      y = runif(max_N_value, 0, max_coord)
    )
    
    points_sf <- sf::st_as_sf(points_df, coords=c("x", "y"), crs=albers_crs, remove=FALSE)
    
    # Extract covariates
    pts_covs <- terra::extract(r_final, points_sf)
    points_df$cell_cov1 <- pts_covs[,2]
    points_df$cell_cov1[is.na(points_df$cell_cov1)] <- 0 # Edge boundary cleanup
    
    ##########
    # 8. ClustGeo Execution
    ##########
    
    # Scale for distance matrices
    scaled_pts <- points_df
    scaled_pts$x <- scale(scaled_pts$x)
    scaled_pts$y <- scale(scaled_pts$y)
    
    geo_dist <- dist(scaled_pts[, c("x", "y")])
    env_dist <- dist(scale(points_df$cell_cov1))
    
    tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = alpha_clustgeo)
    
    K <- round(max_N_value * (kappa_for_clustgeo / 100))
    points_df$site <- as.character(cutree(tree, K))
    
    ##########
    # 9. Build Site Geometries & W Matrix
    ##########
    cat("  Building buffered geometries and computing W...\n")
    
    # Repopulate sf object with site labels
    points_sf$site <- points_df$site
    
    # Voronoi buffering
    site_geoms <- voronoi_clipped_buffers(points_sf, buffer_dist = buffer_m)
    
    # Disjoint Geometry handling
    # Rename x, y to longitude, latitude specifically for disjoint_site_geometries matching logic
    pts_disjoint_input <- points_df
    pts_disjoint_input$longitude <- pts_disjoint_input$x
    pts_disjoint_input$latitude <- pts_disjoint_input$y
    
    split_res <- disjoint_site_geometries(site_geoms, pts_disjoint_input, crs_points = albers_crs)
    final_geoms <- split_res$geoms
    clust_df <- split_res$data 
    
    # Revert to standard x/y for clean tracking
    clust_df$x <- clust_df$longitude; clust_df$y <- clust_df$latitude
    
    # Overlap Matrix
    w_matrix <- generate_overlap_matrix(final_geoms, r_final)
    
    ##########
    # 10. Simulate Site and Observation True States
    ##########
    
    lambda_tilde_i <- as.numeric(w_matrix %*% full_lambda_j)
    
    site_states <- data.frame(
      site = rownames(w_matrix),
      lambda_i = lambda_tilde_i
    )
    site_states$psi_i <- 1 - exp(-site_states$lambda_i)
    site_states$Z_i <- rbinom(nrow(site_states), 1, site_states$psi_i)
    
    # Map Z_i to points
    clust_df <- merge(clust_df, site_states[,c("site", "Z_i")], by="site", all.x=TRUE)
    
    # Generate obs covariates and detections
    clust_df$obs_cov1 <- rnorm(nrow(clust_df))
    logit_p <- true_alphas[1] + true_alphas[2] * clust_df$obs_cov1
    clust_df$p <- plogis(logit_p)
    clust_df$species_observed <- rbinom(nrow(clust_df), 1, clust_df$p * clust_df$Z_i)
    
    # Ensure standard ordering
    clust_df <- clust_df[order(clust_df$checklist_id), ]
    
    ##########
    # 11. INNER LOOP: Ablation over M selected sites
    ##########
    
    available_sites <- unique(clust_df$site)
    
    # Prepare Plotting Lists
    if (sim == 1) {
      plots_cov   <- list()
      plots_abund <- list()
      plots_occ   <- list()
      
      # Base cell dataframe for raster plotting
      cell_df <- as.data.frame(r_final, xy=TRUE)
      names(cell_df)[3] <- "covariate"
    }
    
    for (M_i in M_values_to_test) {
      
      if (M_i > length(available_sites)) {
         cat(sprintf("    Skipping M=%d (Only %d sites generated)\n", M_i, length(available_sites)))
         next
      }
      
      # Randomly Sample M Sites
      selected_sites <- sample(available_sites, M_i, replace = FALSE)
      
      # Subset Data
      sub_clust_df <- clust_df[clust_df$site %in% selected_sites, ]
      sub_w <- w_matrix[selected_sites, , drop = FALSE]
      sub_geoms <- final_geoms[final_geoms$site %in% selected_sites, ]
      
      # Inject into occuPPM formatter
      umf <- prepare_occuPPM_data(
        train_data = sub_clust_df, 
        clustering_df = sub_clust_df[, c("checklist_id", "site")], 
        w_matrix = sub_w, 
        obs_cov_names = c("obs_cov1"), 
        cell_covs = full_cellCovs
      )
      
      # Fit occuPPM
      best_fm <- NULL
      min_nll <- Inf
      n_params <- length(true_alphas) + length(true_betas)
      
      for (rep in 1:n_reps) {
        rand_starts <- runif(n_params, -2, 2) 
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
      
      est_val <- if(is.null(best_fm)) c(NA,NA,NA,NA) else c(coef(best_fm, 'det'), coef(best_fm, 'state'))
      
      # Record Results
      loop_results <- data.frame(
        Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
        True_Value = c(true_alphas, true_betas),
        Estimated_Value = est_val,
        M = M_i, 
        sim_id = sim,
        SAC_Level = sac_level
      )
      results_storage[[length(results_storage) + 1]] <- loop_results
      
      # --- PLOTTING LOGIC (Sim 1 Only) ---
      if (sim == 1) {
        
        # Merge site intensity for visual mapping
        sub_geoms_plot <- merge(sub_geoms, site_states, by="site")
        
        tight_theme <- theme_minimal() + 
          theme(axis.title = element_blank(), plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt"))
        
        # 1. Covariate map with points and borders
        p_cov <- ggplot() +
          geom_raster(data=cell_df, aes(x=x, y=y, fill=covariate)) +
          scale_fill_viridis_c() +
          geom_sf(data=sub_geoms_plot, fill=NA, color="red", linewidth=0.3) +
          geom_point(data=sub_clust_df, aes(x=x, y=y), color="black", size=0.5) +
          coord_sf(expand=FALSE) +
          labs(title=sprintf("Covariate (M=%d)", M_i), fill="Covariate") +
          tight_theme
        
        # 2. Site-Level Abundance (Painted over faded cov layer)
        p_abund <- ggplot() +
          geom_raster(data=cell_df, aes(x=x, y=y, fill=covariate), alpha=0.3, show.legend=FALSE) +
          scale_fill_viridis_c() + 
          ggnewscale::new_scale_fill() +
          geom_sf(data=sub_geoms_plot, aes(fill=lambda_i), color="red", alpha=0.85) +
          scale_fill_viridis_c(option = "magma", name="Abundance") +
          coord_sf(expand=FALSE) +
          labs(title=sprintf("Abundance (M=%d)", M_i)) +
          tight_theme
        
        # 3. Site-Level Occupancy (Painted over faded cov layer)
        p_occ <- ggplot() +
          geom_raster(data=cell_df, aes(x=x, y=y, fill=covariate), alpha=0.3, show.legend=FALSE) +
          scale_fill_viridis_c() + 
          ggnewscale::new_scale_fill() +
          geom_sf(data=sub_geoms_plot, aes(fill=as.factor(Z_i)), color="red", alpha=0.85) +
          scale_fill_manual(values=c("0"="navyblue", "1"="yellow"), name="Occupancy") +
          coord_sf(expand=FALSE) +
          labs(title=sprintf("Occupancy (M=%d)", M_i)) +
          tight_theme
          
        plots_cov[[length(plots_cov)+1]]     <- p_cov
        plots_abund[[length(plots_abund)+1]] <- p_abund
        plots_occ[[length(plots_occ)+1]]     <- p_occ
      }
      
    } # End M Loop
    
    # --- SAVE PLOTS ---
    if (sim == 1) {
      cat(sprintf("  Saving plots for SAC=%s...\n", sac_level))
      
      # Group the plots into columns first without adding themes
      col_cov <- patchwork::wrap_plots(plots_cov, ncol = 1)
      col_abund <- patchwork::wrap_plots(plots_abund, ncol = 1)
      col_occ <- patchwork::wrap_plots(plots_occ, ncol = 1)
      
      # Combine, collect guides, and use '+' instead of '&' to theme the layout
      final_comb_plot <- (col_cov | col_abund | col_occ) + 
        patchwork::plot_layout(guides = "collect") + 
        theme(legend.position = "bottom", legend.direction = "horizontal")
      
      fname <- sprintf("plot_SAC=%s.png", sac_level)
      ggsave(file.path(output_dir, fname), plot=final_comb_plot, dpi=300, width=11, height=20)
    }
    
    gc() # Clean memory
  } # End Sim Loop
} # End SAC Loop


##########
# 12. Save Results & Generate Error Plots
##########

cat("\n--- All Simulations Complete. Saving Results... ---\n")

all_results_df <- do.call(rbind, results_storage)
all_results_df <- all_results_df[!sapply(all_results_df$Parameter, is.null), ]

write.csv(all_results_df, file.path(output_dir, "params_clustgeo.csv"), row.names = FALSE)

# Generate Error Boxplots
all_results_df$Error <- all_results_df$True_Value - all_results_df$Estimated_Value
all_results_df$M_factor <- as.factor(all_results_df$M)
all_results_df$SAC_Level <- factor(all_results_df$SAC_Level, levels = c("Low", "Medium", "High"))

create_error_plot <- function(param_name, title) {
  ggplot(all_results_df[all_results_df$Parameter == param_name, ], 
         aes(x = M_factor, y = Error, fill = SAC_Level)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Low" = "yellow", "Medium" = "orange", "High" = "red")) +
    labs(title = title, x = "M (Sites)", y = "Error (True - Estimate)", fill = "Spatial Autocorrelation") +
    theme_bw() + theme(legend.position = "none")
}

p1 <- create_error_plot("beta (state_int)", "State Intercept")
p2 <- create_error_plot("beta (state_cov1)", "State Slope")
p3 <- create_error_plot("alpha (det_int)", "Observation Intercept")
p4 <- create_error_plot("alpha (det_cov1)", "Observation Slope")

combined_error_plot <- (p1 | p2) / (p3 | p4) +
  patchwork::plot_layout(guides = "collect") + 
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "error_boxplots.png"), plot = combined_error_plot, dpi = 300, width = 10, height = 10)

cat("\n--- Script Finished Successfully ---\n")