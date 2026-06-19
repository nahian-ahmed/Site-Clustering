################################################################
# Simulation Experiments using ClustGeo for Site Generation
#################################################################

###
# 1. SETUP
###

install_now = TRUE
if (install_now){
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  if (!requireNamespace("devtools", quietly = FALSE)) install.packages("devtools")
  suppressMessages(devtools::install_github("anonymous97331/unmarked", ref = "main", force = TRUE))
  if (!requireNamespace("ClustGeo", quietly = FALSE)) install.packages("ClustGeo")
  if (!requireNamespace("sf", quietly = FALSE)) install.packages("sf")
}

library(unmarked)
library(ggplot2)
library(patchwork)
library(terra) 
library(Matrix) 
library(scales)
library(ClustGeo)
library(sf)

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- ClustGeo Parameters ---
kappa_for_clustgeo <- 10  # Percentage of cells to form initial clusters (10% of 40k = 4000)
alpha_for_clustgeo <- 0.50 # Weight between spatial and environmental distance

# --- Simulation repetitions ---
n_sims <- 100 
n_reps <- 30 

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200 
full_n_cells <- full_grid_dim * full_grid_dim # 40000

# --- Observation parameters ---
J_obs <- 3 

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
decay_scale <- 30^2

output_dir <- file.path("output", "simulation_experiments", "clustgeo")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- Simulation Starting (ClustGeo) ---\n")
cat(sprintf("Running %d full simulations per SAC level.\n", n_sims))
cat(sprintf("SAC Levels: %s\n", paste(sac_levels, collapse=", ")))
cat(sprintf("ClustGeo Settings: Kappa = %d%%, Alpha = %.2f\n", kappa_for_clustgeo, alpha_for_clustgeo))
cat(sprintf("Nested Uniform M Sample Sizes: %s\n\n", paste(M_values_to_test, collapse=", ")))

##########
# 3. Base Grid Setup
##########

# Pre-generate coordinates (No physical units, strictly 0 to 200 indices)
full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1

##########
# 4. Pre-generate Fixed Seeds (For Skew)
##########
cat("Pre-generating fixed COVARIATE skew seeds...\n")
cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(n_centers, 0, full_grid_dim), 
    y = runif(n_centers, 0, full_grid_dim)
  )
}

##########
# 5. Initialize Storage
##########
results_storage <- list()

##########
# 6. OUTER LOOPS: SAC -> SIM
##########

for (sac_level in sac_levels) {
  
  current_sigma <- sac_sigmas[[sac_level]]
  
  for (sim in 1:n_sims) {
    
    cat(sprintf("\n=== SAC: %s | Sim %d of %d ===\n", sac_level, sim, n_sims))
    
    ##########
    # 7. Generate Landscape Covariates
    ##########
    
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
    
    full_cellCovs <- data.frame(cell_cov1 = terra::values(r_final, mat=FALSE))
    if(any(is.na(full_cellCovs$cell_cov1))) full_cellCovs$cell_cov1[is.na(full_cellCovs$cell_cov1)] <- 0
    
    ##########
    # 8. ClustGeo Site Generation
    ##########
    cat("  Running ClustGeo on landscape cells...\n")
    
    # Calculate Distances
    # Note: dist() on 40k points uses ~6GB RAM per matrix. 
    env_dist <- dist(scale(full_cellCovs$cell_cov1))
    geo_dist <- dist(cbind(full_cell_col, full_cell_row))
    
    # Run HclustGeo
    tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = alpha_for_clustgeo)
    
    # Cut tree
    init_K <- round(full_n_cells * (kappa_for_clustgeo / 100))
    cluster_assignments <- cutree(tree, init_K)
    
    # Clean up massive dist matrices to free RAM
    rm(env_dist, geo_dist, tree)
    gc()
    
    # --- Split Spatially Disjoint Subclusters ---
    cat("  Splitting disjoint geometries...\n")
    r_clusters <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim)
    terra::values(r_clusters) <- cluster_assignments
    
    # Convert clusters to polygons and separate disconnected parts
    polys <- terra::as.polygons(r_clusters, dissolve=TRUE)
    sf_polys <- sf::st_as_sf(polys) %>% 
      sf::st_cast("MULTIPOLYGON") %>% 
      sf::st_cast("POLYGON", warn = FALSE)
    
    sf_polys$site_id <- 1:nrow(sf_polys)
    M_max <- nrow(sf_polys)
    
    cat(sprintf("  Generated %d spatially contiguous sites (M_max).\n", M_max))
    
    # --- Create W Matrix ---
    r_cells <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim)
    terra::values(r_cells) <- 1:full_n_cells
    
    extracted <- terra::extract(r_cells, terra::vect(sf_polys), cells=FALSE)
    
    full_w <- Matrix::sparseMatrix(
      i = extracted$ID,
      j = extracted[, 2], 
      x = 1,
      dims = c(M_max, full_n_cells)
    )
    
    ##########
    # 9. Generate Biological Data & Observations
    ##########
    
    full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
    full_lambda_j <- exp(full_X_cell %*% true_betas)
    
    full_lambda_tilde_i <- as.numeric(full_w %*% full_lambda_j)
    full_psi_i <- 1 - exp(-full_lambda_tilde_i)
    full_Z_i <- rbinom(M_max, 1, full_psi_i)
    
    full_obs_cov1 <- matrix(rnorm(M_max * J_obs), M_max, J_obs)
    full_obsCovs <- list(obs_cov1 = full_obs_cov1)
    
    full_y <- matrix(NA, M_max, J_obs)
    for (i in 1:M_max) {
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
    # 10. Nested Uniform Sampling Loop
    ##########
    
    # Generate the nested sampling sequence ONCE per simulation
    permuted_sites <- sample(1:M_max)
    
    if (sim == 1) {
      plots_cov   <- list()
      plots_abund <- list()
      plots_occ   <- list()
    }
    
    for (M_i in M_values_to_test) {
      
      if (M_i > M_max) {
        cat(sprintf("  Skipping M=%d (Not enough generated sites)\n", M_i))
        next
      }
      
      # Nested Sample: Always takes the first M_i elements
      selected_site_indices <- permuted_sites[1:M_i]
      
      w_sub <- full_w[selected_site_indices, , drop=FALSE]
      y_sub <- full_y[selected_site_indices, , drop=FALSE]
      obs_cov1_sub <- full_obs_cov1[selected_site_indices, , drop=FALSE]
      obsCovs_sub <- list(obs_cov1 = obs_cov1_sub)
      
      umf <- unmarkedFrameOccuPPM(
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
      
      loop_results <- data.frame(
        Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
        True_Value = c(true_alphas, true_betas),
        Estimated_Value = est_val,
        M = M_i, 
        sim_id = sim,
        SAC_Level = sac_level
      )
      
      results_storage[[length(results_storage) + 1]] <- loop_results
      
      # --- PLOTTING (Sim 1 Only) ---
      if (sim == 1) {
        cell_df <- data.frame(
          x = full_cell_col,
          y = full_cell_row,
          covariate = full_cellCovs$cell_cov1
        )
        
        # We need a cell-to-site map for the *full* generated M_max sites to plot global truth
        cell_df$site_latent_abundance <- 0
        cell_df$site_true_occupancy <- 0
        for(i in 1:M_max) {
           cells_in_i <- which(full_w[i, ] == 1)
           cell_df$site_latent_abundance[cells_in_i] <- full_lambda_tilde_i[i]
           cell_df$site_true_occupancy[cells_in_i] <- full_Z_i[i]
        }
        cell_df$site_true_occupancy <- as.factor(cell_df$site_true_occupancy)
        
        # Get Polygon boundaries for selected sites only
        selected_polys <- sf_polys[sf_polys$site_id %in% selected_site_indices, ]
        
        tight_theme <- theme_minimal() + 
          theme(
            axis.title = element_blank(),
            plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt")
          )
        
        p_cov <- ggplot(cell_df, aes(x=x, y=y, fill=covariate)) +
          geom_raster() +
          scale_fill_viridis_c() +
          coord_fixed(expand=FALSE) +
          geom_sf(data=selected_polys, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          labs(title=sprintf("Covariate (M=%d)", M_i), fill="Covariate") +
          tight_theme
        
        p_abund <- ggplot(cell_df, aes(x=x, y=y, fill=site_latent_abundance)) +
          geom_raster() +
          scale_fill_viridis_c(option = "magma") +
          coord_fixed(expand=FALSE) +
          geom_sf(data=selected_polys, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          labs(title=sprintf("Abundance (M=%d)", M_i), fill="Abundance") +
          tight_theme
        
        p_occ <- ggplot(cell_df, aes(x=x, y=y, fill=site_true_occupancy)) +
          geom_raster() +
          scale_fill_manual(values=c("0"="navyblue", "1"="yellow")) +
          coord_fixed(expand=FALSE) +
          geom_sf(data=selected_polys, color="red", fill=NA, linewidth=0.3, inherit.aes=FALSE) +
          labs(title=sprintf("Occupancy (M=%d)", M_i), fill="Occupancy") +
          tight_theme
          
        plots_cov[[length(plots_cov)+1]]     <- p_cov
        plots_abund[[length(plots_abund)+1]] <- p_abund
        plots_occ[[length(plots_occ)+1]]     <- p_occ
      }
      
    } # End M Loop
    
    # --- SAVE PLOTS (Sim 1 Only) ---
    if (sim == 1) {
      cat(sprintf("\nSaving plots for SAC=%s...\n", sac_level))
      
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
      
      fname <- sprintf("plot_SAC=%s.png", sac_level)
      ggsave(file.path(output_dir, fname), plot=final_comb_plot, dpi=300, width=11, height=20)
    }
    
    gc()
  } # End Sim Loop
} # End SAC Loop


##########
# 11. Save Results & Generate Error Plots 
##########

cat("\n--- All Simulations Complete. Generating Summary Error Boxplots... ---\n")

results_df <- do.call(rbind, results_storage)
results_df <- results_df[!sapply(results_df$Parameter, is.null), ]

write.csv(results_df, file.path(output_dir, "params_clustgeo.csv"), row.names = FALSE)

results_df$Error <- results_df$True_Value - results_df$Estimated_Value
results_df$M_factor <- as.factor(results_df$M)
results_df$SAC_Level <- factor(results_df$SAC_Level, levels = c("Low", "Medium", "High"))

create_error_plot <- function(param_name, title) {
  ggplot(results_df[results_df$Parameter == param_name, ], 
         aes(x = M_factor, y = Error, fill = SAC_Level)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Low" = "yellow", "Medium" = "orange", "High" = "red")) +
    labs(title = title, x = "M (Sites)", y = "Error (True - Estimate)", fill = "Spatial Autocorrelation") +
    theme_bw() + theme(legend.position = "none") # Hide internal legends to collect them cleanly at the end
}

p1 <- create_error_plot("beta (state_int)", "State Intercept")
p2 <- create_error_plot("beta (state_cov1)", "State Slope")
p3 <- create_error_plot("alpha (det_int)", "Observation Intercept")
p4 <- create_error_plot("alpha (det_cov1)", "Observation Slope")

combined_error_plot <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave(file.path(output_dir, "error_boxplots.png"), plot = combined_error_plot, dpi = 300, width = 10, height = 10)

cat(sprintf("Saved summary results and plots to: %s\n", output_dir))
cat("\n--- Script Finished Successfully ---\n")