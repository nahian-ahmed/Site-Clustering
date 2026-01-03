# -----------------------------------------------------------------
# Simulation for occuN model
# Fully simulated experiments with varying Spatial Autocorrelation (SAC)
# AND Sampling Strategies (Random, Positive, Negative, Nonrandom)
# AND Landscape Skew Patterns (Uniform, TopRight, Centered)
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

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Simulation repetitions ---
n_sims <- 100 # Number of full datasets to generate per SAC/Skew level
n_sims <- 3 # For testing

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
M_values_to_test <- c(100, 225, 400, 900)

# --- Sampling Strategies ---
sampling_strategies <- c("Random", "Positive", "Negative", "Nonrandom")

# --- Spatial Autocorrelation (SAC) Settings ---
sac_levels <- c("Low", "Medium", "High")
sac_sigmas <- c(Low = 0, Medium = 1, High = 3) 

# --- NEW: Landscape Skew Patterns ---
# Uniform: Standard SAC (No gradient)
# TopRight: Gradient from Bottom-Left (Low) to Top-Right (High)
# Centered: High values clustered around fixed seed locations ("mountains")
skew_levels <- c("Uniform", "TopRight", "Centered")

# --- Cluster Settings for "Nonrandom" Strategy ---
n_clusters <- 5
cluster_sigma <- 5 # Controls the spread/size of the sampling clusters

cat("--- Simulation Starting ---\n")
cat(sprintf("Running %d full simulations per SAC/Skew level.\n", n_sims))
cat(sprintf("SAC Levels: %s\n", paste(sac_levels, collapse=", ")))
cat(sprintf("Skew Patterns: %s\n", paste(skew_levels, collapse=", ")))
cat(sprintf("Sampling Strategies: %s\n", paste(sampling_strategies, collapse=", ")))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", 
    length(skew_levels) * length(sac_levels) * n_sims * length(M_values_to_test) * length(sampling_strategies) * n_reps))


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

# --- Pre-calculate Site Centroids (for Nonrandom distance calc) ---
site_indices <- 1:full_M
s_rows <- (site_indices - 1) %/% full_n_sites_x + 1
s_cols <- (site_indices - 1) %% full_n_sites_x + 1
site_coords <- data.frame(
  site_id = site_indices,
  x = (s_cols - 1) * site_dim + (site_dim / 2),
  y = (s_rows - 1) * site_dim + (site_dim / 2)
)

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

# 4a. Sampling Cluster Seeds (for "Nonrandom" Strategy)
cat("Pre-generating fixed SAMPLING cluster seeds...\n")
sim_cluster_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  sim_cluster_seeds[[i]] <- list(
    x = runif(n_clusters, 0, full_grid_dim),
    y = runif(n_clusters, 0, full_grid_dim)
  )
}

# 4b. Skew Center Seeds (for "Centered" Skew Pattern)
# These act as the "mountains" for the covariate
cat("Pre-generating fixed COVARIATE skew seeds...\n")
cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(3, 0, full_grid_dim), # 3 mountains
    y = runif(3, 0, full_grid_dim)
  )
}

##########
# 5. Initialize Loop & Storage
##########

total_iterations <- length(skew_levels) * length(sac_levels) * n_sims * length(M_values_to_test) * length(sampling_strategies)
results_list <- vector("list", total_iterations)
results_counter <- 1

output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n--- Starting Main Simulation Loop ---\n")

for (skew in skew_levels) {
    
    cat(sprintf("\n*******************************************************\n"))
    cat(sprintf("*** SKEW PATTERN: %s ***", skew))
    cat(sprintf("\n*******************************************************\n"))
    
    for (sac_level in sac_levels) {
        
        current_sigma <- sac_sigmas[[sac_level]]
        
        for (sim in 1:n_sims) {
            
            cat(sprintf("\n=== Skew: %s | SAC: %s | Sim %d of %d ===\n", skew, sac_level, sim, n_sims))
        
            ##########
            # 6. Create FULL Landscape
            ##########
            
            # --- 6a. Generate Base Spatial Noise (SAC) ---
            r_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                                  xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim,
                                  vals=rnorm(full_n_cells))
            
            if (current_sigma > 0) {
                fw <- terra::focalMat(r_noise, current_sigma, type = "Gauss")
                r_smooth <- terra::focal(r_noise, w = fw, fun = sum, na.rm = TRUE)
            } else {
                r_smooth <- r_noise
            }
            # Normalize noise to Standard Normal (approx) before adding trend
            terra::values(r_smooth) <- as.vector(scale(terra::values(r_smooth)))
            
            # --- 6b. Apply Skew / Trend ---
            r_trend <- r_smooth 
            
            if (skew == "Uniform") {
                # No change, just the GRF
                r_trend <- r_smooth 
                
            } else if (skew == "TopRight") {
                # Create coordinate rasters
                r_x <- terra::init(r_smooth, "x")
                r_y <- terra::init(r_smooth, "y")
                
                # Gradient: Low at (0,0), High at (Max, Max)
                r_gradient <- (r_x + r_y) 
                
                # Normalize gradient to 0-1 range
                v <- terra::values(r_gradient)
                v <- (v - min(v)) / (max(v) - min(v))
                
                # Scale gradient to be strong (-2 to 2) and add to noise
                terra::values(r_gradient) <- v * 4 - 2 
                r_trend <- r_smooth + r_gradient
                
            } else if (skew == "Centered") {
                # Get seeds for this sim
                seeds <- cov_center_seeds[[sim]]
                
                # Initialize trend raster with zeros
                r_centers <- terra::rast(r_smooth)
                terra::values(r_centers) <- 0
                
                # Add Gaussian mountains
                rows <- terra::init(r_smooth, "y")
                cols <- terra::init(r_smooth, "x")
                
                for(k in 1:length(seeds$x)) {
                    # Distance from seed
                    d2 <- (cols - seeds$x[k])^2 + (rows - seeds$y[k])^2
                    # Gaussian decay (broad mountains, sigma=30)
                    r_centers <- r_centers + exp(-d2 / (2 * 30^2))
                }
                
                # Scale centers to be impactful (0 to 3)
                r_trend <- r_smooth + (r_centers * 3)
            }
            
            # --- 6c. Final Standardization ---
            # Important: Ensure mean=0, sd=1 so betas remain comparable across skews
            terra::values(r_trend) <- as.vector(scale(terra::values(r_trend)))
            r_final <- r_trend
            
            full_cellCovs <- data.frame(cell_cov1 = terra::values(r_final, mat=FALSE))
            if(any(is.na(full_cellCovs$cell_cov1))) full_cellCovs$cell_cov1[is.na(full_cellCovs$cell_cov1)] <- 0
            
            # True Lambda
            full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
            full_lambda_j <- exp(full_X_cell %*% true_betas)
            
            # True Site States (for full 1600 sites)
            full_lambda_tilde_i <- as.numeric(full_w %*% full_lambda_j)
            full_psi_i <- 1 - exp(-full_lambda_tilde_i)
            full_Z_i <- rbinom(full_M, 1, full_psi_i)
            
            # --- Calculate Weights for Sampling Strategies ---
            
            # 1. Covariate Based Weights (Positive & Negative)
            site_cov_sums <- as.numeric(full_w %*% full_cellCovs$cell_cov1)
            site_cov_means <- site_cov_sums / (site_dim^2)
            
            # 2. Cluster Based Weights (Nonrandom)
            current_seeds <- sim_cluster_seeds[[sim]]
            
            site_dists <- rep(Inf, full_M)
            for(k in 1:n_clusters){
               d <- sqrt((site_coords$x - current_seeds$x[k])^2 + (site_coords$y - current_seeds$y[k])^2)
               site_dists <- pmin(site_dists, d)
            }
            cluster_weights <- exp(-site_dists^2 / (2 * cluster_sigma^2))


            # --- LOOP OVER SAMPLING STRATEGIES ---
            for (sampling_strat in sampling_strategies) {
                
                if(sim == 1) {
                    plots_cov   <- list()
                    plots_abund <- list()
                    plots_occ   <- list()
                }
                
                cat(sprintf("  >> Strategy: %s\n", sampling_strat))

                # --- SELECT SAMPLING WEIGHTS ---
                if (sampling_strat == "Random") {
                    curr_weights <- NULL
                } else if (sampling_strat == "Positive") {
                    curr_weights <- exp(site_cov_means) 
                } else if (sampling_strat == "Negative") {
                    curr_weights <- exp(-site_cov_means)
                } else if (sampling_strat == "Nonrandom") {
                    curr_weights <- cluster_weights
                }

                # --- LOOP OVER M ---
                for (M_i in M_values_to_test) {
                    
                    # --- 6.1. Select Sites ---
                    if (is.null(curr_weights)) {
                        selected_site_indices <- sample(1:full_M, M_i, replace = FALSE)
                    } else {
                        selected_site_indices <- sample(1:full_M, M_i, replace = FALSE, prob = curr_weights)
                    }
                    
                    selected_site_indices <- sort(selected_site_indices)
                    M <- length(selected_site_indices)
                    
                    # --- 6.2. Subset Data ---
                    w_sub <- full_w[selected_site_indices, , drop=FALSE]
                    Z_sub <- full_Z_i[selected_site_indices]
                    
                    ##########
                    # 8. Simulate Observation Data (y)
                    ##########
                    
                    obs_cov1 <- matrix(rnorm(M * J_obs), M, J_obs)
                    obsCovs <- list(obs_cov1 = obs_cov1)
                    
                    y <- matrix(NA, M, J_obs)
                    for (i in 1:M) {
                        if (Z_sub[i] == 0) {
                            y[i, ] <- 0
                            next
                        }
                        for (k in 1:J_obs) {
                            logit_p_ik <- true_alphas[1] * 1 + true_alphas[2] * obsCovs$obs_cov1[i, k]
                            p_ik <- plogis(logit_p_ik)
                            y[i, k] <- rbinom(1, 1, p_ik)
                        }
                    }
                    
                    ##########
                    # 9. Bundle Data
                    ##########
                    
                    umf <- unmarkedFrameOccuN(
                        y = y,
                        obsCovs = obsCovs,
                        cellCovs = full_cellCovs, 
                        w = w_sub 
                    )
                    
                    ##########
                    # 10. Fit the occuN Model
                    ##########
                    
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
                    
                    # Store Results
                    est_val <- if(is.null(fm)) c(NA,NA,NA,NA) else c(coef(fm, 'det'), coef(fm, 'state'))
                    
                    loop_results <- data.frame(
                        Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
                        True_Value = c(true_alphas, true_betas),
                        Estimated_Value = est_val,
                        M = M, 
                        sim_id = sim,
                        SAC_Level = sac_level,
                        Skew = skew,    # Added Skew
                        Sampling = sampling_strat
                    )
                    
                    results_list[[results_counter]] <- loop_results
                    results_counter <- results_counter + 1
                    
                    ##########
                    # 12. Plotting (Sim 1 Only)
                    ##########
                    if (sim == 1) {
                        cell_df <- data.frame(
                            x = full_cell_col,
                            y = full_cell_row,
                            covariate = full_cellCovs$cell_cov1
                        )
                        
                        cell_df$site_latent_abundance <- full_lambda_tilde_i[full_site_id_for_cell]
                        cell_df$site_true_occupancy <- as.factor(full_Z_i[full_site_id_for_cell])
                        
                        boxes_list <- lapply(selected_site_indices, function(sid) {
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
                            labs(title=sprintf("Covariate (M=%d)", M), fill="Covariate") +
                            tight_theme

                        p_abund <- ggplot(cell_df, aes(x=x, y=y, fill=site_latent_abundance)) +
                            geom_raster() +
                            scale_fill_viridis_c(option = "magma") +
                            coord_fixed(expand=FALSE) +
                            geom_rect(data=site_boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                                      color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
                            labs(title=sprintf("Abundance (M=%d)", M), fill="Abundance") +
                            tight_theme

                        p_occ <- ggplot(cell_df, aes(x=x, y=y, fill=site_true_occupancy)) +
                            geom_raster() +
                            scale_fill_manual(values=c("0"="navyblue", "1"="yellow")) +
                            coord_fixed(expand=FALSE) +
                            geom_rect(data=site_boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                                      color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
                            labs(title=sprintf("Occupancy (M=%d)", M), fill="Occupancy") +
                            tight_theme
                            
                        plots_cov[[length(plots_cov)+1]]     <- p_cov
                        plots_abund[[length(plots_abund)+1]] <- p_abund
                        plots_occ[[length(plots_occ)+1]]     <- p_occ
                    }
                    
                } # End M Loop
                
                # Save Strategy Specific Plots for this SAC/Skew level
                if (sim == 1) {
                    cat(sprintf("Saving plots for Skew=%s, SAC=%s, Sampling=%s...\n", skew, sac_level, sampling_strat))
                    
                    col_cov <- patchwork::wrap_plots(plots_cov, ncol = 1) + 
                      patchwork::plot_layout(guides = "collect") & 
                      theme(legend.position = "bottom", legend.direction = "horizontal")
                      
                    col_abund <- patchwork::wrap_plots(plots_abund, ncol = 1) + 
                      patchwork::plot_layout(guides = "collect") & 
                      theme(legend.position = "bottom", legend.direction = "horizontal")
                    
                    col_occ <- patchwork::wrap_plots(plots_occ, ncol = 1) + 
                      patchwork::plot_layout(guides = "collect") & 
                      theme(legend.position = "bottom", legend.direction = "horizontal")
                    
                    final_comb_plot <- col_cov | col_abund | col_occ
                    
                    # Updated Filename to include Skew
                    fname <- sprintf("plot_Skew=%s_SAC=%s_sampling=%s.png", skew, sac_level, sampling_strat)
                    ggsave(file.path(output_dir, fname), plot=final_comb_plot, dpi=300, width=13, height=18)
                }
                
            } # End Sampling Loop
            
            gc()
            
        } # End Sim Loop
        
    } # End SAC Loop
    
} # End Skew Loop


##########
# 14. Save Aggregate Results & Error Boxplots
##########

cat("\n--- Simulation Study Complete ---\n")

all_results_df <- do.call(rbind, results_list)
all_results_df <- all_results_df[!sapply(all_results_df$Parameter, is.null), ]

write.csv(all_results_df, file.path(output_dir, "params.csv"), row.names = FALSE)

# --- Boxplots Separate by Sampling Strategy ---
# Now iterating over Skew as well to avoid mixing different landscape patterns

for (skew in skew_levels) {
    for (strat in sampling_strategies) {
        
        strat_df <- all_results_df[all_results_df$Sampling == strat & all_results_df$Skew == skew, ]
        if(nrow(strat_df) == 0) next
        
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
        
        fname <- sprintf("error_boxplots_Skew=%s_sampling=%s.png", skew, strat)
        ggsave(file.path(output_dir, fname), plot = combined_error_plot, dpi = 300, width = 10, height = 10)
        cat(sprintf("Saved %s\n", fname))
    }
}

cat("--- Script Finished ---\n")