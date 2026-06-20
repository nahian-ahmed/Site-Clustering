################################################################
# Simulation Experiments with ClustGeo Spatial Clustering
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
library(ClustGeo)
library(dplyr)

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Nested Site Selection ---
M_values_to_test <- c(100, 200, 400, 800, 1600)

max_M <- max(M_values_to_test)

aggreg_factor <- 2
target_K <- max_M

split_factor <- 1

# --- Simulation repetitions ---
n_sims <- 100

# --- Model fitting repetitions ---
n_reps <- 30 

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200 
full_n_cells <- full_grid_dim * full_grid_dim # 40000

# --- ClustGeo Specific Parameters ---
kappa_for_clustgeo <- (max_M / ((full_grid_dim * full_grid_dim)/(aggreg_factor * aggreg_factor)))*100
# kappa_for_clustgeo <- 1  # Used as percentage: 10% of total cells
alpha_for_clustgeo <- 0.9

# --- Observation parameters ---
J_obs <- 3

# --- True parameter values ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
true_betas <- c(beta_int = -5.0, beta_cov = 1.0)

# --- Model settings ---
selected_optimizer <- "nlminb"
PARAM_LOWER <- -20
PARAM_UPPER <- 20
INIT_LOWER <- -5
INIT_UPPER <- 5

# --- Spatial Autocorrelation (SAC) Settings ---
# Fixed to a single "High" value
current_sigma <- 6

# --- Skew Patterns ---
skew <- "Centers"
n_centers <- 1
centers_scale <- 5
decay_scale <- 30^2

output_dir <- file.path("output", "simulation_experiments", "clustgeo")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- Simulation Starting ---\n")
cat(sprintf("Running %d full simulations.\n", n_sims))
cat(sprintf("ClustGeo Settings: Kappa = %d%%, Alpha = %.2f\n", kappa_for_clustgeo, alpha_for_clustgeo))
cat(sprintf("SAC Level (Sigma): %d\n", current_sigma))
cat(sprintf("Nested M Values: %s\n", paste(M_values_to_test, collapse=", ")))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", 
            n_sims * length(M_values_to_test) * n_reps))


##########
# 3. Pre-generate Fixed Seeds & Coordinates
##########

cat("Pre-generating fixed COVARIATE skew seeds...\n")
cov_center_seeds <- vector("list", n_sims)
for(i in 1:n_sims){
  cov_center_seeds[[i]] <- list(
    x = runif(n_centers, 0, full_grid_dim), 
    y = runif(n_centers, 0, full_grid_dim)
  )
}

# Base Abstract Raster (No physical units)
r_base <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                      xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim)

# Cell Coordinates (centers)
full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1
cell_centers <- terra::xyFromCell(r_base, 1:full_n_cells)
centers_sf <- sf::st_as_sf(as.data.frame(cell_centers), coords=c("x", "y"))
sf::st_crs(centers_sf) <- NA

##########
# 4. Initialize Storage
##########

results_storage <- list()

##########
# 5. OUTER LOOP: SIMULATIONS
##########

for (sim in 1:n_sims) {
  
  cat(sprintf("\n=== Sim %d of %d ===\n", sim, n_sims))
  
  ##########
  # 6. Generate Landscape Data
  ##########
  
  # --- 6a. Generate Base Spatial Noise (SAC) ---
  r_noise <- terra::rast(r_base)
  terra::values(r_noise) <- rnorm(full_n_cells)
  
  if (current_sigma > 0) {
    fw <- terra::focalMat(r_noise, current_sigma, type = "Gauss")
    r_smooth <- terra::focal(r_noise, w = fw, fun = sum, na.rm = TRUE)
  } else {
    r_smooth <- r_noise
  }
  terra::values(r_smooth) <- as.vector(scale(terra::values(r_smooth)))
  
  # --- 6b. Apply Skew / Trend (Centers) ---
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
  # 7. ClustGeo Workflow (2x2 Aggregation -> Cluster -> Disaggregate)
  ##########
  
  cat("  Running ClustGeo spatial clustering...\n")
  
  # 1. Aggregate 4 cells (2x2)
  r_agg <- terra::aggregate(r_final, fact=aggreg_factor, fun=mean) 
  df_agg <- as.data.frame(r_agg, xy=TRUE)
  names(df_agg) <- c("x", "y", "cov")
  
  # 2. ClustGeo execution
  env_dist <- dist(df_agg$cov)
  geo_dist <- dist(df_agg[, c("x", "y")])
  tree <- ClustGeo::hclustgeo(env_dist, geo_dist, alpha = alpha_for_clustgeo)
  
  # K parameter: user specified kappa=10 => 10% of 40,000 cells = 4,000 clusters
  # K_target <- (kappa_for_clustgeo / 100) * full_n_cells
  K_target <- as.integer(target_K /split_factor)
  cluster_ids <- cutree(tree, k = K_target)
  
  # 3. Map back and Disaggregate
  r_clust_agg <- terra::rast(r_agg)
  terra::values(r_clust_agg) <- cluster_ids
  r_clust <- terra::disagg(r_clust_agg, fact=2, method="near")
  
  # 4. Spatially disjoint splitting using SF
  polys <- terra::as.polygons(r_clust, dissolve=TRUE)
  polys_sf <- sf::st_as_sf(polys)
  sf::st_crs(polys_sf) <- NA
  
  # Split MULTIPOLYGON into POLYGON to separate disconnected subclusters
  polys_split <- suppressWarnings(
      sf::st_cast(sf::st_cast(polys_sf, "MULTIPOLYGON"), "POLYGON")
  )
  
  full_M_max <- nrow(polys_split)
  polys_split$site_id <- 1:full_M_max
  
  cat(sprintf("  ClustGeo yielded %d final disjoint site geometries.\n", full_M_max))
  
  # Map cells back to new disjoint site geometries
  inter <- sf::st_intersects(centers_sf, polys_split)
  full_site_id_for_cell <- sapply(inter, function(x) x[1])
  
  # Create W Matrix
  full_w <- Matrix::sparseMatrix(
    i = full_site_id_for_cell,
    j = 1:full_n_cells,
    x = 1,
    dims = c(full_M_max, full_n_cells)
  )
  
  ##########
  # 8. Generate Biological Data & Observations
  ##########
  
  # True expected abundance per cell
  full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
  full_lambda_j <- exp(full_X_cell %*% true_betas)
  
  # Site-level abundance and occupancy
  full_lambda_tilde_i <- as.numeric(full_w %*% full_lambda_j)
  full_psi_i <- 1 - exp(-full_lambda_tilde_i)
  full_Z_i <- rbinom(full_M_max, 1, full_psi_i)
  
  # Observations (J=3)
  full_obs_cov1 <- matrix(rnorm(full_M_max * J_obs), full_M_max, J_obs)
  full_obsCovs <- list(obs_cov1 = full_obs_cov1)
  
  full_y <- matrix(NA, full_M_max, J_obs)
  for (i in 1:full_M_max) {
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
  # 9. Nested Sampling & Model Fitting
  ##########
  
  # Create a random permutation for perfectly nested sampling
  site_permutation <- sample(1:full_M_max, full_M_max, replace = FALSE)
  
  if (sim == 1) {
    all_plots_list <- list()
  }
  
  for (M_i in M_values_to_test) {
    
    # Select perfectly nested subset
    selected_site_indices <- sort(site_permutation[1:M_i])
    M <- length(selected_site_indices)
    
    # Subset Data
    w_sub <- full_w[selected_site_indices, , drop=FALSE]
    y_sub <- full_y[selected_site_indices, , drop=FALSE]
    obs_cov1_sub <- full_obs_cov1[selected_site_indices, , drop=FALSE]
    obsCovs_sub <- list(obs_cov1 = obs_cov1_sub)
    
    # Fit Model
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
      sim_id = sim
    )
    
    results_storage[[length(results_storage) + 1]] <- loop_results
    
    # --- PLOTTING LOGIC (Only for sim == 1) ---
    if (sim == 1) {
      
      cell_df <- data.frame(
        x = full_cell_col,
        y = full_cell_row,
        covariate = full_cellCovs$cell_cov1
      )
      cell_df$site_latent_abundance <- full_lambda_tilde_i[full_site_id_for_cell]
      cell_df$site_true_occupancy <- as.factor(full_Z_i[full_site_id_for_cell])
      
      selected_polys <- polys_split[polys_split$site_id %in% selected_site_indices, ]
      
      tight_theme <- theme_minimal() + 
        theme(
          axis.title = element_blank(),
          plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt")
        )
      
      # Note: Using coord_sf(expand=FALSE, datum=NA) avoids coordinate system warnings
      p_cov <- ggplot(cell_df, aes(x=x, y=y, fill=covariate)) +
        geom_raster() +
        scale_fill_viridis_c() +
        geom_sf(data=selected_polys, fill=NA, color="red", linewidth=0.5, inherit.aes=FALSE) +
        coord_sf(expand=FALSE, datum=NA) +
        labs(title=sprintf("Covariate (M=%d)", M), fill="Covariate") +
        tight_theme
      
      p_abund <- ggplot(cell_df, aes(x=x, y=y, fill=site_latent_abundance)) +
        geom_raster() +
        scale_fill_viridis_c(option = "magma") +
        geom_sf(data=selected_polys, fill=NA, color="red", linewidth=0.5, inherit.aes=FALSE) +
        coord_sf(expand=FALSE, datum=NA) +
        labs(title=sprintf("Abundance (M=%d)", M), fill="Abundance") +
        tight_theme
      
      p_occ <- ggplot(cell_df, aes(x=x, y=y, fill=site_true_occupancy)) +
        geom_raster() +
        scale_fill_manual(values=c("0"="navyblue", "1"="yellow")) +
        geom_sf(data=selected_polys, fill=NA, color="red", linewidth=0.5, inherit.aes=FALSE) +
        coord_sf(expand=FALSE, datum=NA) +
        labs(title=sprintf("Occupancy (M=%d)", M), fill="Occupancy") +
        tight_theme
        
      # --- NEW LOGIC: Format legends flat and ONLY on the last row ---
      if (M_i == max(M_values_to_test)) {
        # For the last row: position bottom, format horizontally, put title on top
        p_cov <- p_cov + 
          theme(legend.position = "bottom") +
          guides(fill = guide_colorbar(direction = "horizontal", barwidth = 10, title.position = "top", title.hjust = 0.5))
        
        p_abund <- p_abund + 
          theme(legend.position = "bottom") +
          guides(fill = guide_colorbar(direction = "horizontal", barwidth = 10, title.position = "top", title.hjust = 0.5))
        
        p_occ <- p_occ + 
          theme(legend.position = "bottom", legend.direction = "horizontal") +
          guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
      } else {
        # For all other rows: remove the legend completely
        p_cov <- p_cov + theme(legend.position = "none")
        p_abund <- p_abund + theme(legend.position = "none")
        p_occ <- p_occ + theme(legend.position = "none")
      }
      
      all_plots_list <- c(all_plots_list, list(p_cov, p_abund, p_occ))
    }
    
  } # End M Loop
  
  # --- SAVE PLOTS (After M loop, for Sim 1) ---
  if (sim == 1) {
    cat("\nSaving landscape plots...\n")
    
    # Just wrap the plots. The bottom three already have horizontal legends attached!
    combined_plot <- patchwork::wrap_plots(all_plots_list, 
                                           nrow = length(M_values_to_test), 
                                           ncol = 3)
                                               
    fname <- "plot.png"
    ggsave(file.path(output_dir, fname), plot=combined_plot, dpi=300, width=9, height=16)
  }
  
  gc()
  
} # End Sim Loop


##########
# 10. Save Results & Generate Error Plots
##########

cat("\n--- All Simulations Complete. Saving Results... ---\n")

all_results_df <- do.call(rbind, results_storage)

if (!is.null(all_results_df) && nrow(all_results_df) > 0) {
    
    all_results_df <- all_results_df[!sapply(all_results_df$Parameter, is.null), ]
    
    # Save CSV
    write.csv(all_results_df, file.path(output_dir, "params_clustgeo.csv"), row.names = FALSE)
    
    # Boxplots
    all_results_df$Error <- all_results_df$True_Value - all_results_df$Estimated_Value
    all_results_df$M_factor <- as.factor(all_results_df$M)
    
    create_error_plot <- function(param_name, title) {
      ggplot(all_results_df[all_results_df$Parameter == param_name, ], 
             aes(x = M_factor, y = Error)) +
        # --- CHANGED: color="red" for the outline, fill=NA for no fill ---
        geom_boxplot(color = "black", fill = "grey", outlier.size = 0.5, alpha = 0.25) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
        labs(title = title, x = "M (Sites)", y = "Error (True - Estimate)") +
        theme_bw()
    }
    
    p1 <- create_error_plot("beta (state_int)", "State Intercept")
    p2 <- create_error_plot("beta (state_cov1)", "State Slope")
    p3 <- create_error_plot("alpha (det_int)", "Observation Intercept")
    p4 <- create_error_plot("alpha (det_cov1)", "Observation Slope")
    
    combined_error_plot <- (p1 | p2) / (p3 | p4)
    
    ggsave(file.path(output_dir, "error_boxplots.png"), plot = combined_error_plot, dpi = 300, width = 8, height = 8)
    
    cat("Saved results and error boxplots.\n")
}

cat("\n--- Script Finished Successfully ---\n")