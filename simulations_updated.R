################################################################
# Updated Simulation Experiments: Sampling Extents
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

##########
# 2. Set Simulation Parameters
##########

set.seed(123) 

# --- Simulation repetitions ---
n_sims <- 3
n_reps <- 30 

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200 
full_n_cells <- full_grid_dim * full_grid_dim # 40000

# --- Observation parameters ---
J_obs <- 3 

# --- True parameter values ---
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)
# true_betas <- c(beta_int = -5.0, beta_cov = 1.0)
true_betas <- c(beta_int = -5.0, beta_cov = 1.0)

# --- Model settings ---
selected_optimizer <- "nlminb"
PARAM_LOWER <- -20
PARAM_UPPER <- 20

# --- Fixed SAC and Skew ---
sac_sigma <- 3 # Medium SAC
n_centers <- 1
centers_scale <- 5
decay_scale <- 30^2

# --- 3 Scenarios: Sampling Extents ---
extents <- c("Small" = 1600, "Medium" = 400, "Large" = 100)

# FORCE absolute path
output_dir <- file.path(getwd(), "output", "simulation_experiments", "updated")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- Simulation Starting ---\n")
cat(sprintf("Saving all outputs to exactly:\n -> %s\n\n", normalizePath(output_dir)))
cat(sprintf("Running %d simulations for 3 Sampling Extents.\n", n_sims))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", n_sims * length(extents) * n_reps))


##########
# 3. Pre-generate Coordinates & Seeds
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
# 4. Main Simulation Loop
##########

all_results <- list()
plot_data <- list()

for (sim in 1:n_sims) {
  cat(sprintf("\n=== Sim %d of %d ===\n", sim, n_sims))
  
  # --- 4a. Base Covariate ---
  r_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                         xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim,
                         vals=rnorm(full_n_cells))
  
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
  
  # --- 4b. True Cell States (Realized N and 0/1) ---
  full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
  full_lambda_j <- exp(full_X_cell %*% true_betas)
  
  full_N_j <- rpois(full_n_cells, full_lambda_j)
  full_Z_j <- factor(ifelse(full_N_j > 0, 1, 0), levels = c("0", "1"))
  
  if (sim == 1) {
    plot_data[["Cell"]] <- data.frame(
      x = full_cell_col, y = full_cell_row,
      covariate = full_cellCovs$cell_cov1,
      abundance = full_N_j,
      occupancy = full_Z_j
    )
  }
  
  # --- 4c. Generate Covariate-Biased Site Geometries ---
  if (sim == 1) cat("Generating Density-Weighted Voronoi Site Geometries...\n")
  site_definitions <- list()
  
  # Create a highly peaked probability weight based on the covariate
  # Scale by 2.5 to create strong hotspots
  # prob_weights <- exp(full_cellCovs$cell_cov1 * 2.5) 
  prob_weights <- exp(full_cellCovs$cell_cov1 * 0.5)
  
  for (ext_name in names(extents)) {
    K <- extents[[ext_name]]
    
    # 1. Sample K seed cells strongly biased towards high covariate values
    seed_idx <- sample(1:full_n_cells, K, prob = prob_weights)
    seed_pts <- cell_coords[seed_idx, ]
    
    # 2. Fast Base-R Voronoi Assignment (Assigns all 40,000 cells to nearest seed)
    site_ids <- rep(1, full_n_cells)
    min_dists <- rep(Inf, full_n_cells)
    for (k in 1:K) {
      dists_sq <- (cell_coords$x - seed_pts$x[k])^2 + (cell_coords$y - seed_pts$y[k])^2
      update_idx <- dists_sq < min_dists
      site_ids[update_idx] <- k
      min_dists[update_idx] <- dists_sq[update_idx]
    }
    
    w <- Matrix::sparseMatrix(
      i = site_ids,
      j = 1:full_n_cells,
      x = 1,
      dims = c(K, full_n_cells)
    )
    
    # Extract Polygons ONLY on sim 1 to save massive compute time on sims 2+
    site_sf <- NULL
    if (sim == 1) {
      xyz_df <- data.frame(x = full_cell_col, y = full_cell_row, z = site_ids)
      r_sites <- terra::rast(xyz_df, type = "xyz")
      site_polys <- terra::as.polygons(r_sites, dissolve = TRUE)
      site_sf <- sf::st_as_sf(site_polys)
      colnames(site_sf)[1] <- "site"
    }
    
    site_definitions[[ext_name]] <- list(
      K = K,
      site_ids = site_ids,
      w = w,
      site_sf = site_sf
    )
  }
  
  # --- 4d. Iterate Over Extents ---
  for (ext_name in names(extents)) {
    cat(sprintf("  - Extent: %s\n", ext_name))
    
    def <- site_definitions[[ext_name]]
    M <- def$K
    w <- def$w
    site_ids <- def$site_ids
    
    # Calculate Site-level True States
    lambda_tilde_i <- as.numeric(w %*% full_lambda_j)
    N_i <- rpois(M, lambda_tilde_i)
    Z_i <- factor(ifelse(N_i > 0, 1, 0), levels = c("0", "1"))
    
    # Generate Observations (y)
    obs_cov1 <- matrix(rnorm(M * J_obs), M, J_obs)
    obsCovs <- list(obs_cov1 = obs_cov1)
    y <- matrix(NA, M, J_obs)
    
    for (i in 1:M) {
      if (Z_i[i] == "0") {
        y[i, ] <- 0
        next
      }
      for (k in 1:J_obs) {
        logit_p_ik <- true_alphas[1] * 1 + true_alphas[2] * obsCovs$obs_cov1[i, k]
        y[i, k] <- rbinom(1, 1, plogis(logit_p_ik))
      }
    }
    
    if (sim == 1) {
      agg_cov <- as.numeric((w %*% full_cellCovs$cell_cov1) / rowSums(w))
      
      sf_data <- def$site_sf
      sf_data <- sf_data[order(sf_data$site), ] 
      sf_data$covariate <- agg_cov
      sf_data$abundance <- N_i
      sf_data$occupancy <- Z_i
      
      plot_data[[ext_name]] <- sf_data
    }
    
    # Fit occuPPM Model
    umf <- unmarkedFrameOccuPPM(
      y = y, obsCovs = obsCovs, cellCovs = full_cellCovs, w = w 
    )
    
    best_fm <- NULL
    min_nll <- Inf
    n_params <- length(true_alphas) + length(true_betas)
    
    for (rep in 1:n_reps) {
      rand_starts <- runif(n_params, -2, 2) 
      fm_rep <- try(occuPPM(
        formula = ~obs_cov1 ~ cell_cov1,
        data = umf, starts = rand_starts, se = TRUE,
        method = selected_optimizer, lower = PARAM_LOWER, upper = PARAM_UPPER 
      ), silent = TRUE)
      
      if (inherits(fm_rep, "try-error")) next
      curr_nll <- fm_rep@negLogLike
      if (curr_nll < min_nll) {
        min_nll <- curr_nll
        best_fm <- fm_rep
      }
    } 
    
    est_val <- if(is.null(best_fm)) c(NA,NA,NA,NA) else c(coef(best_fm, 'det'), coef(best_fm, 'state'))
    
    all_results[[length(all_results) + 1]] <- data.frame(
      Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
      True_Value = c(true_alphas, true_betas),
      Estimated_Value = est_val,
      Extent = ext_name,
      sim_id = sim
    )
  }
}

##########
# 5. Generate and Save 4x3 Plot (Sampling Extents)
##########
cat("\nGenerating 4x3 Spatial Plot (sampling_extents.png)...\n")

cov_limits <- range(c(plot_data$Cell$covariate, plot_data$Small$covariate, plot_data$Medium$covariate, plot_data$Large$covariate), na.rm=TRUE)
abund_limits <- range(c(plot_data$Cell$abundance, plot_data$Small$abundance, plot_data$Medium$abundance, plot_data$Large$abundance), na.rm=TRUE)

ggplot2::theme_set(ggplot2::theme_minimal() + ggplot2::theme(
  legend.position = "bottom",
  legend.justification = "center",
  legend.box = "horizontal",
  legend.box.just = "center",
  legend.spacing.x = ggplot2::unit(2.5, "cm"), 
  legend.box.margin = ggplot2::margin(t = 5, l = 15), 
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
  direction = "horizontal",
  title.position = "top",
  title.hjust = 0.5,
  barwidth = ggplot2::unit(5.0, "cm"),
  barheight = ggplot2::unit(0.6, "cm")
)

guide_disc <- ggplot2::guide_legend(
  direction = "horizontal",
  title.position = "top",
  title.hjust = 0.5,
  keywidth = ggplot2::unit(1.0, "cm"), 
  keyheight = ggplot2::unit(0.6, "cm")
)

build_row <- function(data_name, row_title, show_titles=FALSE) {
  df <- plot_data[[data_name]]
  
  if (data_name == "Cell") {
    p1 <- ggplot(df, aes(x=x, y=y, fill=covariate)) + geom_raster() +
      scale_fill_viridis_c(name="Covariate", limits=cov_limits, guide=guide_cont) + base_theme +
      labs(y = row_title, x = NULL) + coord_fixed(expand=FALSE)
      
    p2 <- ggplot(df, aes(x=x, y=y, fill=abundance)) + geom_raster() +
      scale_fill_viridis_c(option="magma", name="Abundance", limits=abund_limits, guide=guide_cont) + base_theme +
      labs(y = NULL, x = NULL) + coord_fixed(expand=FALSE)
      
    p3 <- ggplot(df, aes(x=x, y=y, fill=occupancy)) + geom_raster() +
      scale_fill_manual(values=c("0"="#440154FF", "1"="#FDE725FF"), name="           Occupancy           ", drop=FALSE, guide=guide_disc) + base_theme +
      labs(y = NULL, x = NULL) + coord_fixed(expand=FALSE)
      
  } else {
    p1 <- ggplot(df) + geom_sf(aes(fill=covariate), color="black", linewidth=0.1) +
      scale_fill_viridis_c(name="Covariate", limits=cov_limits, guide=guide_cont) + base_theme +
      labs(y = row_title, x = NULL) + coord_sf(expand=FALSE)
      
    p2 <- ggplot(df) + geom_sf(aes(fill=abundance), color="black", linewidth=0.1) +
      scale_fill_viridis_c(option="magma", name="Abundance", limits=abund_limits, guide=guide_cont) + base_theme +
      labs(y = NULL, x = NULL) + coord_sf(expand=FALSE)
      
    p3 <- ggplot(df) + geom_sf(aes(fill=occupancy), color="black", linewidth=0.1, show.legend=FALSE) +
      scale_fill_manual(values=c("0"="#440154FF", "1"="#FDE725FF"), name="           Occupancy           ", drop=FALSE, guide=guide_disc) + base_theme +
      labs(y = NULL, x = NULL) + coord_sf(expand=FALSE)
  }
  
  if(show_titles) {
    p1 <- p1 + ggtitle("Covariate")
    p2 <- p2 + ggtitle("Abundance")
    p3 <- p3 + ggtitle("Occupancy")
  }
  
  return(list(p1, p2, p3))
}

row1 <- build_row("Cell", "Simulated Species", show_titles=TRUE)
row2 <- build_row("Small", "Small Sampling Extents\n(M = 1600)")
row3 <- build_row("Medium", "Medium Sampling Extents\n(M = 400)")
row4 <- build_row("Large", "Large Sampling Extents\n(M = 100)")

comb_plot <- patchwork::wrap_plots(c(row1, row2, row3, row4), ncol=3) + 
  patchwork::plot_layout(guides="collect")

ggsave(file.path(output_dir, "sampling_extents.png"), plot=comb_plot, width=10, height=14, dpi=300)

##########
# 6. Process Results & Generate Error Boxplots
##########
cat("Saving Results & Generating Error Boxplots...\n")

res_df <- do.call(rbind, all_results)
res_df <- res_df[!is.na(res_df$Estimated_Value), ]

res_df$Error <- res_df$True_Value - res_df$Estimated_Value
res_df$Extent <- factor(res_df$Extent, levels = c("Small", "Medium", "Large"))

write.csv(res_df, file.path(output_dir, "params_updated.csv"), row.names = FALSE)

create_error_plot <- function(param_name, title) {
  ggplot(res_df[res_df$Parameter == param_name, ], 
         aes(x = Extent, y = Error, fill = Extent)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Small" = "lightblue", "Medium" = "steelblue", "Large" = "navy")) +
    labs(title = title, x = "Sampling Extent", y = "Error (True - Estimate)") +
    theme_bw() + theme(legend.position = "none")
}

p_beta0 <- create_error_plot("beta (state_int)", "State Intercept")
p_beta1 <- create_error_plot("beta (state_cov1)", "State Slope")
p_alpha0 <- create_error_plot("alpha (det_int)", "Observation Intercept")
p_alpha1 <- create_error_plot("alpha (det_cov1)", "Observation Slope")

combined_error_plot <- (p_beta0 | p_beta1) / (p_alpha0 | p_alpha1)
ggsave(file.path(output_dir, "error_boxplots.png"), plot = combined_error_plot, dpi = 300, width = 9, height = 9)

cat("--- Script Finished Successfully ---\n")