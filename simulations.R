# -----------------------------------------------------------------
# Simulation for occuN model
# Fully simulated experiments with varying Spatial Autocorrelation (SAC)
# Fixed Skew Pattern: Centers
# Variable Sampling Strategies: Uniform, Positive (High Cov), Negative (Low Cov)
# -----------------------------------------------------------------

###
# 1. SETUP
###

install_now = FALSE
if (install_now){
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  if (!requireNamespace("devtools", quietly = FALSE)) install.packages("devtools")
  if (!requireNamespace("terra", quietly = FALSE)) install.packages("terra")
  suppressMessages(devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE))
}

library(unmarked)
library(ggplot2)
library(patchwork)
library(terra) 
library(Matrix) 
library(scales) 

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Simulation repetitions ---
n_sims <- 100 
# n_sims <- 3 # FOR DEBUGGING

# --- Model fitting repetitions ---
n_reps <- 30 

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200 
full_n_cells <- full_grid_dim * full_grid_dim # 40000

# --- Site parameters ---
site_dim <- 5 # Sites are 5x5 cell blocks
full_n_sites_x <- full_grid_dim / site_dim # 40
full_n_sites_y <- full_grid_dim / site_dim # 40
full_M <- full_n_sites_x * full_n_sites_y # 1600

# --- Observation parameters ---
J_obs <- 3 

# --- True parameter values ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
true_betas <- c(beta_int = -5.0, beta_cov = 1.0)


# --- Model settings ---
selected_optimizer <- "nlminb"
# Optimization bounds
PARAM_LOWER <- -20
PARAM_UPPER <- 20

# --- Ablation Study Parameters ---
M_values_to_test <- c(100, 200, 400, 800, 1600)

# --- Sampling Strategies ---
sampling_strategies <- c("Uniform", "Positive", "Negative")

# --- Weighted Sampling Parameters ---
n_sampling_clusters <- 1   # Number of hotspots to select
sampling_sd <- 5           # SD for Gaussian smoothing (in site units, ~25 cells)

# --- Spatial Autocorrelation (SAC) Settings ---
sac_levels <- c("Low", "Medium", "High") 
# sac_sigmas <- c(Low = 0, Medium = 5, High = 15)
sac_sigmas <- c(Low = 0, Medium = 10, High = 30) 

# --- Skew Patterns ---
# Fixed to Centers
skew <- "Centers"

# --- Centers for "Centers" Skew ---
n_centers <- 1
centers_scale <- 5
decay_scale <- 30^2

output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- Simulation Starting ---\n")
cat(sprintf("Running %d full simulations per SAC level.\n", n_sims))
cat(sprintf("SAC Levels: %s\n", paste(sac_levels, collapse=", ")))
cat(sprintf("Strategies (Inner Loop): %s\n", paste(sampling_strategies, collapse=", ")))
cat(sprintf("Weighted Sampling: %d clusters, SD=%d\n", n_sampling_clusters, sampling_sd))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", 
            length(sampling_strategies) * length(sac_levels) * n_sims * length(M_values_to_test) * n_reps))


##########
# 3. Define FULL Sites & Create Weight Matrix (w)
##########

cat("Generating static full landscape geometry (1600 sites from 40000 cells)...\n")

# Assign each cell to a site
full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1
full_site_row <- (full_cell_row - 1) %/% site_dim + 1
full_site_col <- (full_cell_col - 1) %/% site_dim + 1

full_site_id_for_cell <- (full_site_row - 1) * full_n_sites_x + full_site_col

# Full weight matrix w (1600 x 40000)
full_w <- Matrix::sparseMatrix(
  i = full_site_id_for_cell,
  j = 1:full_n_cells,
  x = 1,
  dims = c(full_M, full_n_cells)
)

# Site coordinates (1 to 40)
site_coords_x <- rep(1:full_n_sites_x, times = full_n_sites_y)
site_coords_y <- rep(1:full_n_sites_y, each = full_n_sites_x)

# Helper to get site boxes for plotting
get_site_box <- function(site_id, site_dim, n_sites_x) {
  s_row <- (site_id - 1) %/% n_sites_x + 1
  s_col <- (site_id - 1) %% n_sites_x + 1
  
  xmin <- (s_col - 1) * site_dim + 0.5
  xmax <- xmin + site_dim
  ymin <- (s_row - 1) * site_dim + 0.5
  ymax <- ymin + site_dim
  return(c(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax))
}

##########
# 4. Pre-generate Fixed Seeds
##########

# 4b. Skew Center Seeds (for "Centers" Skew Pattern)
cat("Pre-generating fixed COVARIATE skew seeds...\n")
cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(n_centers, 0, full_grid_dim), # 3 mountains
    y = runif(n_centers, 0, full_grid_dim)
  )
}

##########
# 5. Initialize Storage
##########

# We will store results in a list of lists, keyed by Strategy Name
# e.g. results_storage[["Uniform"]] will hold all rows for uniform sampling
results_storage <- list()
for(strat in sampling_strategies) {
  results_storage[[strat]] <- list()
}

##########
# 6. OUTER LOOPS: SAC -> SIM
##########

for (sac_level in sac_levels) {
  
  current_sigma <- sac_sigmas[[sac_level]]
  
  for (sim in 1:n_sims) {
    
    cat(sprintf("\n=== SAC: %s | Sim %d of %d ===\n", sac_level, sim, n_sims))
    
    ##########
    # 7. Generate Landscape & Biological Data (Once per Sim)
    ##########
    
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
    
    # --- 7b. Apply Skew / Trend (Centers) ---
    seeds <- cov_center_seeds[[sim]] # Get fixed seeds for this sim
    
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
    
    full_cellCovs <- data.frame(cell_cov1 = terra::values(r_final, mat=FALSE))
    if(any(is.na(full_cellCovs$cell_cov1))) full_cellCovs$cell_cov1[is.na(full_cellCovs$cell_cov1)] <- 0
    
    # --- 7c. Generate True States (Z) ---
    full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
    full_lambda_j <- exp(full_X_cell %*% true_betas)
    
    full_lambda_tilde_i <- as.numeric(full_w %*% full_lambda_j)
    full_psi_i <- 1 - exp(-full_lambda_tilde_i)
    full_Z_i <- rbinom(full_M, 1, full_psi_i)
    
    # --- 7d. Generate Observations (y) ---
    full_obs_cov1 <- matrix(rnorm(full_M * J_obs), full_M, J_obs)
    full_obsCovs <- list(obs_cov1 = full_obs_cov1)
    
    full_y <- matrix(NA, full_M, J_obs)
    for (i in 1:full_M) {
      if (full_Z_i[i] == 0) {
        full_y[i, ] <- 0
        next
      }
      for (k in 1:J_obs) {
        logit_p_ik <- true_alphas[1] * 1 + true_alphas[2] * full_obsCovs$obs_cov1[i, k]
        p_ik <- plogis(logit_p_ik)
        full_y[i, k] <- rbinom(1, 1, p_ik)
      }
    }
    
    ##########
    # 8. INNER LOOP: SAMPLING STRATEGIES
    ##########
    # Slice the *same* landscape/observations using different strategies
    
    for (sampling_strat in sampling_strategies) {
      
      # --- 8a. Calculate Sampling Weights ---
      
      site_sampling_weights <- rep(1, full_M)
      
      if (sampling_strat != "Uniform") {
        
        # 1. Aggregate Covariate to Site Level
        site_cov_vals <- tapply(full_cellCovs$cell_cov1, full_site_id_for_cell, mean)
        
        # 2. Calculate Distance to Skew Centers
        # Convert seed coordinates (grid units 0-200) to site units (0-40)
        # Note: Site coordinates are indices 1..40
        seed_x_site <- seeds$x / site_dim
        seed_y_site <- seeds$y / site_dim
        
        # For simplicity with multiple centers, we take distance to the closest center
        # (Though simulation uses n_centers=1 usually)
        dist_to_nearest_center <- rep(Inf, full_M)
        
        for(k in seq_along(seed_x_site)) {
            d <- sqrt((site_coords_x - seed_x_site[k])^2 + (site_coords_y - seed_y_site[k])^2)
            dist_to_nearest_center <- pmin(dist_to_nearest_center, d)
        }
        
        # 3. Determine Cluster Centers
        cluster_center_ids <- NULL
        
        if (sampling_strat == "Positive") {
          # POSITIVE: Force selection of sites CLOSEST to the Skew Center
          # "On top of each other" -> Deterministic selection based on distance
          
          # Sort sites by distance to skew center
          sorted_indices <- order(dist_to_nearest_center)
          # Pick top N closest sites
          cluster_center_ids <- sorted_indices[1:n_sampling_clusters]
          
        } else if (sampling_strat == "Negative") {
          # NEGATIVE: Probabilistic, favoring Low Covariate & Far Distance
          center_probs <- exp(-site_cov_vals) * (dist_to_nearest_center + 1e-6)
          cluster_center_ids <- sample(1:full_M, size = n_sampling_clusters, prob = center_probs, replace = FALSE)
        }
        
        # 4. Grow Clusters (Gaussian Smoothing)
        prob_surface <- rep(0, full_M)
        
        center_coords_df <- data.frame(
          id = cluster_center_ids,
          x = site_coords_x[cluster_center_ids],
          y = site_coords_y[cluster_center_ids]
        )
        
        for (i in 1:full_M) {
          sx <- site_coords_x[i]
          sy <- site_coords_y[i]
          
          w_sum <- 0
          for (c_idx in 1:nrow(center_coords_df)) {
            cx <- center_coords_df$x[c_idx]
            cy <- center_coords_df$y[c_idx]
            dist_sq <- (sx - cx)^2 + (sy - cy)^2
            w_sum <- w_sum + exp(-dist_sq / (2 * sampling_sd^2))
          }
          prob_surface[i] <- w_sum
        }
        
        site_sampling_weights <- prob_surface
      }
      
      # --- 8b. Generate Plots (Sim 1 Only) ---
      if (sim == 1) {
          # Select sites for M=max just for visualization context, or M=400
          vis_M <- 400
          if (sum(site_sampling_weights > 0) < vis_M) {
               vis_indices <- sample(1:full_M, vis_M, replace = FALSE)
          } else {
               vis_indices <- sample(1:full_M, vis_M, replace = FALSE, prob = site_sampling_weights)
          }
          
          # Prepare Plot Data
          cell_df <- data.frame(
            x = full_cell_col,
            y = full_cell_row,
            covariate = full_cellCovs$cell_cov1
          )
          cell_df$site_latent_abundance <- full_lambda_tilde_i[full_site_id_for_cell]
          cell_df$site_true_occupancy <- as.factor(full_Z_i[full_site_id_for_cell])
          
          boxes_list <- lapply(vis_indices, function(sid) {
            coords <- get_site_box(sid, site_dim, full_n_sites_x)
            data.frame(xmin=coords['xmin'], xmax=coords['xmax'], 
                       ymin=coords['ymin'], ymax=coords['ymax'])
          })
          site_boxes <- do.call(rbind, boxes_list)
          
          tight_theme <- theme_minimal() + 
            theme(
              axis.title = element_blank(),
              plot.margin = margin(t=0, r=0, b=0, l=0, unit="pt")
            )
          
          p_cov <- ggplot(cell_df, aes(x=x, y=y, fill=covariate)) +
            geom_raster() +
            scale_fill_viridis_c() +
            coord_fixed(expand=FALSE) +
            geom_rect(data=site_boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                      color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
            labs(title=sprintf("Covariate (Vis M=%d)", vis_M), fill="Covariate") +
            tight_theme
          
          p_abund <- ggplot(cell_df, aes(x=x, y=y, fill=site_latent_abundance)) +
            geom_raster() +
            scale_fill_viridis_c(option = "magma") +
            coord_fixed(expand=FALSE) +
            geom_rect(data=site_boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                      color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
            labs(title="Abundance", fill="Abundance") +
            tight_theme
          
          p_occ <- ggplot(cell_df, aes(x=x, y=y, fill=site_true_occupancy)) +
            geom_raster() +
            scale_fill_manual(values=c("0"="navyblue", "1"="yellow")) +
            coord_fixed(expand=FALSE) +
            geom_rect(data=site_boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                      color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
            labs(title="Occupancy", fill="Occupancy") +
            tight_theme
          
          final_comb_plot <- p_cov | p_abund | p_occ
          
          # Plot Filename: Include SAC and Strategy
          fname <- sprintf("plot_SAC=%s_sampling=%s.png", sac_level, sampling_strat)
          ggsave(file.path(output_dir, fname), plot=final_comb_plot, dpi=300, width=11, height=6)
      }
      
      # --- 8c. Loop Over M (Sample Size) ---
      for (M_i in M_values_to_test) {
        
        # Select Sites
        if (sum(site_sampling_weights > 0) < M_i) {
             selected_site_indices <- sample(1:full_M, M_i, replace = FALSE)
        } else {
             selected_site_indices <- sample(1:full_M, M_i, replace = FALSE, prob = site_sampling_weights)
        }
        selected_site_indices <- sort(selected_site_indices)
        M <- length(selected_site_indices)
        
        # Subset Data
        w_sub <- full_w[selected_site_indices, , drop=FALSE]
        y_sub <- full_y[selected_site_indices, , drop=FALSE]
        obs_cov1_sub <- full_obs_cov1[selected_site_indices, , drop=FALSE]
        obsCovs_sub <- list(obs_cov1 = obs_cov1_sub)
        
        # Fit Model
        umf <- unmarkedFrameOccuN(
          y = y_sub,
          obsCovs = obsCovs_sub,
          cellCovs = full_cellCovs, 
          w = w_sub 
        )
        
        best_fm <- NULL
        min_nll <- Inf
        n_params <- length(true_alphas) + length(true_betas)
        
        for (rep in 1:n_reps) {
          rand_starts <- runif(n_params, -2, 2) 
          fm_rep <- try(occuN(
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
          M = M, 
          sim_id = sim,
          SAC_Level = sac_level,
          Skew = skew,
          Sampling = sampling_strat
        )
        
        # Append to strategy specific list
        # (This is more efficient than growing a huge dataframe inside the inner loop)
        results_storage[[sampling_strat]][[length(results_storage[[sampling_strat]]) + 1]] <- loop_results
        
      } # End M Loop
      
    } # End Strategy Loop (Inner)
    
    gc()
    
  } # End Sim Loop (Mid)
  
} # End SAC Loop (Outer)


##########
# 9. Save Results & Generate Error Plots (Per Strategy)
##########

cat("\n--- All Simulations Complete. Saving Results per Strategy... ---\n")

for (strat in sampling_strategies) {
  
  # Combine list of dataframes into one
  strat_results_df <- do.call(rbind, results_storage[[strat]])
  
  if (!is.null(strat_results_df) && nrow(strat_results_df) > 0) {
      
      # Clean
      strat_results_df <- strat_results_df[!sapply(strat_results_df$Parameter, is.null), ]
      
      # Save CSV
      csv_name <- sprintf("params_sampling=%s.csv", strat)
      write.csv(strat_results_df, file.path(output_dir, csv_name), row.names = FALSE)
      
      # Boxplots
      strat_df <- strat_results_df
      strat_df$Error <- strat_df$True_Value - strat_df$Estimated_Value
      strat_df$M_factor <- as.factor(strat_df$M)
      strat_df$SAC_Level <- factor(strat_df$SAC_Level, levels = c("Low", "Medium", "High"))
      
      create_error_plot <- function(param_name, title) {
        ggplot(strat_df[strat_df$Parameter == param_name, ], 
               aes(x = M_factor, y = Error, fill = SAC_Level)) +
          geom_boxplot(outlier.size = 0.5) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
          scale_fill_manual(values = c("Low" = "yellow", "Medium" = "orange", "High" = "red")) +
          labs(title = title, x = "M (Sites)", y = "Error", fill = "Spatial Autocorrelation") +
          theme_bw() + theme(legend.position = "none")
      }
      
      p1 <- create_error_plot("beta (state_int)", "State Intercept")
      p2 <- create_error_plot("beta (state_cov1)", "State Slope")
      p3 <- create_error_plot("alpha (det_int)", "Detection Intercept")
      p4 <- create_error_plot("alpha (det_cov1)", "Detection Slope")
      
      combined_error_plot <- (p1 | p2) / (p3 | p4) +
        plot_layout(guides = "collect") & theme(legend.position = "bottom")
      
      fname_plot <- sprintf("error_boxplots_sampling=%s.png", strat)
      ggsave(file.path(output_dir, fname_plot), plot = combined_error_plot, dpi = 300, width = 10, height = 10)
      
      cat(sprintf("Saved results and plots for strategy: %s\n", strat))
  }
}

cat("\n--- Script Finished Successfully ---\n")