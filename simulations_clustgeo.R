################################################################
# Simulation Experiments: "Virtual Observer" & ClustGeo
# Fixes Target Visits per Site & Scales K proportionally
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
library(ClustGeo)
library(dplyr)
library(tidyr)
library(sf)

# Load existing project helpers
source(file.path("R", "utils.R"))
source(file.path("R", "model_helpers.R"))
source(file.path("R", "analysis_helpers.R"))

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Simulation repetitions ---
n_sims <- 3
n_reps <- 30 # Number of random-start repetitions for each model fit

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200 
full_n_cells <- full_grid_dim * full_grid_dim # 40000

# To integrate with real spatial functions, we treat each cell as 1km x 1km
cell_size_m <- 1000 

# --- True parameter values ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
true_betas <- c(beta_int = -5.0, beta_cov = 1.0)

selected_optimizer <- "nlminb"
PARAM_LOWER <- -20
PARAM_UPPER <- 20


PARAM_LOWER <- -10
PARAM_UPPER <- 10
INIT_LOWER <- -2
INIT_UPPER <- 2

# --- "Virtual Observer" Study Parameters ---
N_values_to_test <- c(800, 1200, 1600, 2000, 2400) # Number of observations (checklists)

# We want an average of 10 visits per site to ensure model convergence
target_visits_per_site <- 10 
buffer_dist_m <- 2 * cell_size_m # Buffer of 2 cells

# --- Sampling Strategies ---
sampling_strategies <- c("Uniform", "Positive", "Negative")

# --- Weighted Sampling Parameters ---
n_sampling_clusters <- 1   
sampling_sd_cells <- 25 # SD for Gaussian smoothing in cell units 

# --- Spatial Autocorrelation (SAC) Settings ---
sac_levels <- c("Low", "Medium", "High") 
sac_sigmas <- c(Low = 6, Medium = 12, High = 24) 

# --- Skew Patterns ---
skew <- "Centers"
n_centers <- 1
centers_scale <- 5
decay_scale <- 30^2

output_dir <- file.path("output", "simulation_experiments", "clustgeo")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- ClustGeo Virtual Observer Simulation Starting ---\n")
cat(sprintf("Running %d full simulations per SAC level.\n", n_sims))
cat(sprintf("Observations (N): %s\n", paste(N_values_to_test, collapse=", ")))
cat(sprintf("Target average visits per site: %d\n", target_visits_per_site))

##########
# 3. Define Continuous Landscape Reference
##########

# Coordinates in Cell Units
full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1

# Base Reference Raster (Pseudo-Mercator so we can use meters naturally)
pseudo_crs <- "EPSG:3857"
ref_raster <- terra::rast(
  nrows = full_grid_dim, ncols = full_grid_dim, 
  xmin = 0, xmax = full_grid_dim * cell_size_m, 
  ymin = 0, ymax = full_grid_dim * cell_size_m,
  crs = pseudo_crs
)
terra::values(ref_raster) <- 1:full_n_cells

##########
# 4. Pre-generate Fixed Seeds
##########

cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(n_centers, 0, full_grid_dim), 
    y = runif(n_centers, 0, full_grid_dim)
  )
}

results_storage <- list()
for(strat in sampling_strategies) results_storage[[strat]] <- list()

##########
# 6. OUTER LOOPS: SAC -> SIM
##########

for (sac_level in sac_levels) {
  current_sigma <- sac_sigmas[[sac_level]]
  
  for (sim in 1:n_sims) {
    cat(sprintf("\n=== SAC: %s | Sim %d of %d ===\n", sac_level, sim, n_sims))
    
    # --- 7a. Generate Base Spatial Noise (SAC) ---
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
    
    # --- 7b. Apply Skew / Trend ---
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
    
    # --- 7c. True Latent Cell Abundance ---
    full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
    full_lambda_j <- exp(full_X_cell %*% true_betas)

    ##########
    # 8. INNER LOOP: SAMPLING STRATEGIES
    ##########
    
    for (sampling_strat in sampling_strategies) {
      
      # --- 8a. Calculate Cell Sampling Weights ---
      cell_sampling_weights <- rep(1, full_n_cells)
      
      if (sampling_strat != "Uniform") {
        if (sampling_strat == "Positive") {
          sampling_probs <- exp(full_cellCovs$cell_cov1 * 5)
        } else {
          sampling_probs <- exp(-full_cellCovs$cell_cov1 * 5)
        }
        cluster_center_ids <- sample(1:full_n_cells, size = n_sampling_clusters, prob = sampling_probs, replace = FALSE)
        
        prob_surface <- rep(0, full_n_cells)
        for (i in 1:full_n_cells) {
          sx <- full_cell_col[i]; sy <- full_cell_row[i]
          w_sum <- 0
          for (c_id in cluster_center_ids) {
            cx <- full_cell_col[c_id]; cy <- full_cell_row[c_id]
            dist_sq <- (sx - cx)^2 + (sy - cy)^2
            w_sum <- w_sum + exp(-dist_sq / (2 * sampling_sd_cells^2))
          }
          prob_surface[i] <- w_sum
        }
        cell_sampling_weights <- prob_surface
      }
      
      if (sim == 1) {
        plots_cov <- list(); plots_abund <- list(); plots_occ <- list()
      }
      
      # --- 8b. Loop Over N (Observations) ---
      for (N_i in N_values_to_test) {
        
        # 1. Drop "Virtual Checklists" (Sample with replacement)
        if (sum(cell_sampling_weights > 0) < N_i) {
             sampled_cell_ids <- sample(1:full_n_cells, N_i, replace = TRUE)
        } else {
             sampled_cell_ids <- sample(1:full_n_cells, N_i, replace = TRUE, prob = cell_sampling_weights)
        }
        
        unique_cells <- unique(sampled_cell_ids)
        
        # 2. ClustGeo Clustering
        target_K <- round(N_i / target_visits_per_site)
        K <- max(2, min(length(unique_cells), target_K))
        
        locs_m_x <- (full_cell_col[unique_cells] - 0.5) * cell_size_m
        locs_m_y <- (full_cell_row[unique_cells] - 0.5) * cell_size_m
        
        geo_dist <- dist(data.frame(x = locs_m_x, y = locs_m_y))
        env_dist <- dist(full_cellCovs$cell_cov1[unique_cells])
        
        tree <- hclustgeo(env_dist, geo_dist, alpha = 0.5)
        site_assignments <- as.character(cutree(tree, K))
        
        # 3. Create Dataframes (with jitter to prevent GEOS crash)
        jitter_x <- runif(length(unique_cells), -5, 5)
        jitter_y <- runif(length(unique_cells), -5, 5)
        
        locs_df <- data.frame(
          locality_id = unique_cells,
          longitude = locs_m_x + jitter_x,
          latitude = locs_m_y + jitter_y,
          site = site_assignments
        )
        
        checklists_df <- data.frame(
          checklist_id = 1:N_i,
          locality_id = sampled_cell_ids,
          longitude = locs_df$longitude[match(sampled_cell_ids, locs_df$locality_id)],
          latitude = locs_df$latitude[match(sampled_cell_ids, locs_df$locality_id)],
          site = locs_df$site[match(sampled_cell_ids, locs_df$locality_id)]
        )
        
        # 4. Generate Geometries & Handle Disjoint Sites
        pts_sf <- sf::st_as_sf(locs_df, coords = c("longitude", "latitude"), crs = pseudo_crs)
        
        site_geoms <- voronoi_clipped_buffers(pts_sf, buffer_dist = buffer_dist_m)
        split_res <- disjoint_site_geometries(site_geoms, checklists_df, crs_points = pseudo_crs)
        
        site_geoms <- split_res$geoms
        checklists_df <- split_res$data
        
        Final_M <- nrow(site_geoms)
        
        # 5. Generate W Matrix
        w_matrix <- generate_overlap_matrix(site_geoms, ref_raster)
        
        # 6. Simulate True State (Z_i) from W Matrix & Cell Density
        lambda_tilde <- as.numeric(w_matrix %*% full_lambda_j)
        names(lambda_tilde) <- rownames(w_matrix)
        
        psi_i <- 1 - exp(-lambda_tilde)
        Z_i <- rbinom(Final_M, 1, psi_i)
        names(Z_i) <- rownames(w_matrix)
        
        # 7. Simulate Observations (Y) for Checklists
        checklists_df$Z <- Z_i[as.character(checklists_df$site)]
        checklists_df <- checklists_df[!is.na(checklists_df$Z), ] # Drop orphans
        
        checklists_df$obs_cov1 <- rnorm(nrow(checklists_df))
        checklists_df$p <- plogis(true_alphas[1] + true_alphas[2] * checklists_df$obs_cov1)
        checklists_df$species_observed <- rbinom(nrow(checklists_df), 1, checklists_df$Z * checklists_df$p)
        
        # 8. Prepare occuPPM Data
        umf <- prepare_occuPPM_data(
          train_data = checklists_df, 
          clustering_df = checklists_df, 
          w_matrix = w_matrix, 
          obs_cov_names = c("obs_cov1"), 
          cell_covs = full_cellCovs
        )
        
        # 9. Fit Model
        best_fm <- NULL
        min_nll <- Inf
        n_params <- length(true_alphas) + length(true_betas)
        
        for (rep in 1:n_reps) {
          rand_starts <- runif(n_params, INIT_LOWER, INIT_UPPER) 
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
          
          if (fm_rep@negLogLike < min_nll) {
            min_nll <- fm_rep@negLogLike
            best_fm <- fm_rep
          }
        } 
        
        est_val <- if(is.null(best_fm)) c(NA,NA,NA,NA) else c(coef(best_fm, 'det'), coef(best_fm, 'state'))
        
        # 10. Record Result
        loop_results <- data.frame(
          Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
          True_Value = c(true_alphas, true_betas),
          Estimated_Value = est_val,
          N_obs = N_i, 
          M_clusters = Final_M,
          sim_id = sim,
          SAC_Level = sac_level,
          Skew = skew,
          Sampling = sampling_strat
        )
        results_storage[[sampling_strat]][[length(results_storage[[sampling_strat]]) + 1]] <- loop_results
        
        # 11. Plotting Logic (Sim 1 Only)
        if (sim == 1) {
          
          plot_geoms <- site_geoms
          plot_geoms$site_char <- as.character(plot_geoms$site)
          plot_geoms$Abundance <- lambda_tilde[plot_geoms$site_char]
          plot_geoms$Occupancy <- as.factor(Z_i[plot_geoms$site_char])

          cell_df <- data.frame(
            x = (full_cell_col - 0.5) * cell_size_m,
            y = (full_cell_row - 0.5) * cell_size_m,
            covariate = full_cellCovs$cell_cov1
          )
          
          # FIX: We put the legend settings directly into the ggplot theme here!
          tight_theme <- theme_minimal() + 
            theme(axis.title = element_blank(), axis.text = element_blank(),
                  plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt"),
                  legend.position = "bottom", legend.direction = "horizontal")
          
          p_cov <- ggplot() +
            geom_raster(data=cell_df, aes(x=x, y=y, fill=covariate)) +
            scale_fill_viridis_c() +
            geom_sf(data=plot_geoms, fill=NA, color="red", linewidth=0.3) +
            geom_point(data=checklists_df, aes(x=longitude, y=latitude), size=0.5, color="white", alpha=0.6) +
            coord_sf(expand=FALSE) +
            labs(title=sprintf("Covariate (N=%d)", N_i), fill="Covariate") +
            tight_theme
            
          p_abund <- ggplot() +
            geom_sf(data=plot_geoms, aes(fill=Abundance), color="red", linewidth=0.3) +
            scale_fill_viridis_c(option = "magma") +
            coord_sf(expand=FALSE) +
            labs(title=sprintf("Abundance (N=%d)", N_i), fill="Abundance") +
            tight_theme

          p_occ <- ggplot() +
            geom_sf(data=plot_geoms, aes(fill=Occupancy), color="red", linewidth=0.3) +
            scale_fill_manual(values=c("0"="navyblue", "1"="yellow")) +
            coord_sf(expand=FALSE) +
            labs(title=sprintf("Occupancy (N=%d)", N_i), fill="Occupancy") +
            tight_theme

          plots_cov[[length(plots_cov)+1]] <- p_cov
          plots_abund[[length(plots_abund)+1]] <- p_abund
          plots_occ[[length(plots_occ)+1]] <- p_occ
        }
        
      } # End N Loop
      
      # Save Combined 3-Column Plots
      if (sim == 1) {
        cat(sprintf("\nSaving plots for SAC=%s, Sampling=%s...\n", sac_level, sampling_strat))
        
        col_cov <- patchwork::wrap_plots(plots_cov, ncol = 1) 
        col_abund <- patchwork::wrap_plots(plots_abund, ncol = 1) 
        col_occ <- patchwork::wrap_plots(plots_occ, ncol = 1) 
        
        # FIX: Simply assemble and collect. Patchwork inherits the "bottom" setting automatically!
        final_comb_plot <- (col_cov | col_abund | col_occ) + patchwork::plot_layout(guides = "collect") 
        
        fname <- sprintf("plot_SAC=%s_sampling=%s.png", sac_level, sampling_strat)
        ggsave(file.path(output_dir, fname), plot=final_comb_plot, dpi=300, width=11, height=20)
      }
      
    } # End Strategy Loop
    gc()
  } # End Sim Loop
} # End SAC Loop


##########
# 9. Save Results & Generate Error Plots
##########

cat("\n--- All Simulations Complete. Saving Results... ---\n")

for (strat in sampling_strategies) {
  strat_results_df <- do.call(rbind, results_storage[[strat]])
  
  if (!is.null(strat_results_df) && nrow(strat_results_df) > 0) {
      strat_results_df <- strat_results_df[!sapply(strat_results_df$Parameter, is.null), ]
      write.csv(strat_results_df, file.path(output_dir, sprintf("params_sampling=%s.csv", strat)), row.names = FALSE)
      
      strat_df <- strat_results_df
      strat_df$Error <- strat_df$True_Value - strat_df$Estimated_Value
      
      strat_df$N_factor <- factor(strat_df$N_obs, levels = N_values_to_test)
      strat_df$SAC_Level <- factor(strat_df$SAC_Level, levels = c("Low", "Medium", "High"))
      
      create_error_plot <- function(param_name, title) {
        ggplot(strat_df[strat_df$Parameter == param_name, ], 
               aes(x = N_factor, y = Error, fill = SAC_Level)) +
          geom_boxplot(outlier.size = 0.5) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
          scale_fill_manual(values = c("Low" = "yellow", "Medium" = "orange", "High" = "red")) +
          labs(title = title, x = "N (Observations)", y = "Error (True - Estimate)", fill = "Spatial Autocorrelation") +
          theme_bw() + 
          theme(legend.position = "bottom") # FIX: Base plot gets the legend position!
      }
      
      p1 <- create_error_plot("beta (state_int)", "State Intercept")
      p2 <- create_error_plot("beta (state_cov1)", "State Slope")
      p3 <- create_error_plot("alpha (det_int)", "Observation Intercept")
      p4 <- create_error_plot("alpha (det_cov1)", "Observation Slope")
      
      # FIX: Assemble and collect. No theme operators needed!
      combined_error_plot <- (p1 | p2) / (p3 | p4) + patchwork::plot_layout(guides = "collect") 
      
      ggsave(file.path(output_dir, sprintf("error_boxplots_sampling=%s.png", strat)), plot = combined_error_plot, dpi = 300, width = 10, height = 10)
  }
}

cat("\n--- Script Finished Successfully ---\n")