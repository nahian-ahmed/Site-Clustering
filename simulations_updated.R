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
n_sims <- 5 
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

# --- Fixed SAC and Skew ---
sac_sigma <- 3 # Medium SAC
n_centers <- 1
centers_scale <- 5
decay_scale <- 30^2

# --- 3 Scenarios: Sampling Extents ---
# Defined by the number of sites the landscape is divided into
extents <- c("Small" = 1600, "Medium" = 400, "Large" = 100)

output_dir <- file.path("output", "simulation_experiments", "updated")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("--- Simulation Starting ---\n")
cat(sprintf("Output Directory: %s\n", output_dir))
cat(sprintf("Running %d simulations for 3 Sampling Extents.\n", n_sims))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", n_sims * length(extents) * n_reps))


##########
# 3. Pre-generate Landscape Geometries (Irregular Sites)
##########

cat("Pre-generating irregular, contiguous site geometries using K-means...\n")

full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1
cell_coords <- data.frame(x = full_cell_col, y = full_cell_row)

site_definitions <- list()

for (ext_name in names(extents)) {
  K <- extents[[ext_name]]
  
  # K-means naturally partitions the coordinates into contiguous Voronoi regions
  # This creates non-overlapping, contiguous sites of varying sizes and irregular shapes
  set.seed(42 + K) 
  km <- kmeans(cell_coords, centers = K, iter.max = 100)
  site_ids <- km$cluster
  
  # Build sparse W matrix (Rows = Sites, Cols = Cells)
  w <- Matrix::sparseMatrix(
    i = site_ids,
    j = 1:full_n_cells,
    x = 1,
    dims = c(K, full_n_cells)
  )
  
  # Polygonize sites for plotting boundaries
  r_sites <- terra::rast(nrows = full_grid_dim, ncols = full_grid_dim, 
                         xmin = 0.5, xmax = full_grid_dim + 0.5, 
                         ymin = 0.5, ymax = full_grid_dim + 0.5)
  terra::values(r_sites) <- site_ids
  site_polys <- terra::as.polygons(r_sites, dissolve = TRUE)
  site_sf <- sf::st_as_sf(site_polys)
  
  site_definitions[[ext_name]] <- list(
    K = K,
    site_ids = site_ids,
    w = w,
    site_sf = site_sf
  )
}

# Pre-generate Skew Center Seeds
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
plot_data <- list() # Store data for plotting the 1st simulation

for (sim in 1:n_sims) {
  cat(sprintf("\n=== Sim %d of %d ===\n", sim, n_sims))
  
  # --- 4a. Generate Base Spatial Covariate (SAC = Medium) ---
  r_noise <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                         xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim,
                         vals=rnorm(full_n_cells))
  
  fw <- terra::focalMat(r_noise, sac_sigma, type = "Gauss")
  r_smooth <- terra::focal(r_noise, w = fw, fun = sum, na.rm = TRUE)
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
  
  full_cellCovs <- data.frame(cell_cov1 = terra::values(r_trend, mat=FALSE))
  if(any(is.na(full_cellCovs$cell_cov1))) full_cellCovs$cell_cov1[is.na(full_cellCovs$cell_cov1)] <- 0
  
  # True Cell States (for baseline/plotting)
  full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
  full_lambda_j <- exp(full_X_cell %*% true_betas)
  full_psi_j <- 1 - exp(-full_lambda_j)
  
  # Store 1st row data for Sim 1 Plotting
  if (sim == 1) {
    plot_data[["Cell"]] <- data.frame(
      x = full_cell_col, y = full_cell_row,
      covariate = full_cellCovs$cell_cov1,
      abundance = full_lambda_j,
      occupancy = full_psi_j
    )
  }
  
  # --- 4b. Iterate Over Extents ---
  for (ext_name in names(extents)) {
    cat(sprintf("  - Extent: %s\n", ext_name))
    
    def <- site_definitions[[ext_name]]
    M <- def$K
    w <- def$w
    site_ids <- def$site_ids
    
    # Calculate Site-level True States
    lambda_tilde_i <- as.numeric(w %*% full_lambda_j)
    psi_i <- 1 - exp(-lambda_tilde_i)
    Z_i <- rbinom(M, 1, psi_i)
    
    # Generate Observations (y)
    obs_cov1 <- matrix(rnorm(M * J_obs), M, J_obs)
    obsCovs <- list(obs_cov1 = obs_cov1)
    y <- matrix(NA, M, J_obs)
    
    for (i in 1:M) {
      if (Z_i[i] == 0) {
        y[i, ] <- 0
        next
      }
      for (k in 1:J_obs) {
        logit_p_ik <- true_alphas[1] * 1 + true_alphas[2] * obsCovs$obs_cov1[i, k]
        y[i, k] <- rbinom(1, 1, plogis(logit_p_ik))
      }
    }
    
    # Store aggregated data mapped to cells for Sim 1 Plotting
    if (sim == 1) {
      # Mean covariate per site mapped back to the raster
      agg_cov <- as.numeric((w %*% full_cellCovs$cell_cov1) / rowSums(w))
      
      plot_data[[ext_name]] <- data.frame(
        x = full_cell_col, y = full_cell_row,
        covariate = agg_cov[site_ids],
        abundance = lambda_tilde_i[site_ids],
        occupancy = psi_i[site_ids]
      )
    }
    
    # Fit occuPPM Model
    umf <- unmarkedFrameOccuPPM(
      y = y,
      obsCovs = obsCovs,
      cellCovs = full_cellCovs, 
      w = w 
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

# Moved the legend theme settings directly into base_theme to avoid the '&' error
base_theme <- theme_minimal() + theme(
  axis.text = element_blank(), axis.ticks = element_blank(),
  panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2, "pt"),
  plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
  axis.title.y = element_text(size = 12, face = "bold", angle = 90, vjust = 0.5),
  legend.position = "bottom", 
  legend.direction = "horizontal", 
  legend.key.width = unit(1.5, "cm")
)

build_row <- function(data_name, row_title, show_titles=FALSE, sf_data=NULL) {
  df <- plot_data[[data_name]]
  
  # P1: Covariate
  p1 <- ggplot(df, aes(x=x, y=y, fill=covariate)) + geom_raster() +
    scale_fill_viridis_c(name="Covariate") + base_theme +
    labs(y = row_title, x = NULL)
  if(show_titles) p1 <- p1 + ggtitle("Covariate")
  
  # P2: Abundance
  p2 <- ggplot(df, aes(x=x, y=y, fill=abundance)) + geom_raster() +
    scale_fill_viridis_c(option="magma", name="Abundance") + base_theme +
    labs(y = NULL, x = NULL)
  if(show_titles) p2 <- p2 + ggtitle("Abundance")
  
  # P3: Occupancy
  p3 <- ggplot(df, aes(x=x, y=y, fill=occupancy)) + geom_raster() +
    scale_fill_viridis_c(option="plasma", name="Occupancy") + base_theme +
    labs(y = NULL, x = NULL)
  if(show_titles) p3 <- p3 + ggtitle("Occupancy")
  
  # Conditionally apply SF boundaries and proper coordinate limits 
  # to avoid the "Coordinate system already present" warnings
  if (!is.null(sf_data)) {
    p1 <- p1 + geom_sf(data=sf_data, fill=NA, color="black", linewidth=0.1, inherit.aes=FALSE) + coord_sf(expand=FALSE)
    p2 <- p2 + geom_sf(data=sf_data, fill=NA, color="black", linewidth=0.1, inherit.aes=FALSE) + coord_sf(expand=FALSE)
    p3 <- p3 + geom_sf(data=sf_data, fill=NA, color="black", linewidth=0.1, inherit.aes=FALSE) + coord_sf(expand=FALSE)
  } else {
    p1 <- p1 + coord_fixed(expand=FALSE)
    p2 <- p2 + coord_fixed(expand=FALSE)
    p3 <- p3 + coord_fixed(expand=FALSE)
  }
  
  return(list(p1, p2, p3))
}

row1 <- build_row("Cell", "Simulated Species", show_titles=TRUE)
row2 <- build_row("Small", "Sampling Extent 1\n(Small)", sf_data=site_definitions[["Small"]]$site_sf)
row3 <- build_row("Medium", "Sampling Extent 2\n(Medium)", sf_data=site_definitions[["Medium"]]$site_sf)
row4 <- build_row("Large", "Sampling Extent 3\n(Large)", sf_data=site_definitions[["Large"]]$site_sf)

# Assemble using patchwork
comb_plot <- patchwork::wrap_plots(c(row1, row2, row3, row4), ncol=3) + 
  patchwork::plot_layout(guides="collect")

ggsave(file.path(output_dir, "sampling_extents.png"), plot=comb_plot, width=10, height=13, dpi=300)

##########
# 6. Process Results & Generate Error Boxplots
##########
cat("Saving Results & Generating Error Boxplots...\n")

res_df <- do.call(rbind, all_results)
res_df <- res_df[!is.na(res_df$Estimated_Value), ]

# Calculate Error (True - Estimated)
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