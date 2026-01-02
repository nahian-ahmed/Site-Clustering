
# -----------------------------------------------------------------
# Simulation for occuN model
# Fully simulated experiments with varying Spatial Autocorrelation (SAC)
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
library(terra) # Required for spatial autocorrelation generation

##########
# 2. Set Simulation Parameters
##########

set.seed(123) # For reproducibility

# --- Simulation repetitions ---
n_sims <- 100 # Number of full datasets to generate per SAC level
# n_sims <- 3
# --- Model fitting repetitions ---
n_reps <- 30 # Number of random-start repetitions for each model fit

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200 # Landscape is 200x200 cells
full_n_cells <- full_grid_dim * full_grid_dim # Total number of cells (40000)

# --- Site parameters (reference clustering) ---
site_dim <- 5 # Sites are 5x5 cell blocks (5-cellSq)
full_n_sites_x <- full_grid_dim / site_dim # 40
full_n_sites_y <- full_grid_dim / site_dim # 40
full_M <- full_n_sites_x * full_n_sites_y # Total number of sites (1600)

# --- Observation parameters ---
J_obs <- 3 # Number of surveys per site

# --- True parameter values ---
# Detection (alphas) for formula ~obs_cov1
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)

# State (betas) for formula ~cell_cov1
true_betas <- c(beta_int = -5.0, beta_cov = 1.0) 

# --- Model settings ---
selected_optimizer <- "nlminb"

# --- Ablation Study Parameters ---
M_values_to_test <- c(100, 225, 400, 900, 1600)

# --- Spatial Autocorrelation (SAC) Settings ---
# We test 3 levels of spatial autocorrelation.
# We use a Gaussian filter with varying sigma to induce correlation.
sac_levels <- c("Low", "Medium", "High")
# Sigma values for Gaussian smoothing (0 = random/Low, higher = more smooth)
sac_sigmas <- c(Low = 0, Medium = 5, High = 15) 
sac_sigmas <- c(Low = 0, Medium = 5, High = 10) 


cat("--- Simulation Starting ---\n")
cat(sprintf("Running %d full simulations per SAC level.\n", n_sims))
cat(sprintf("SAC Levels: %s\n", paste(sac_levels, collapse=", ")))
cat(sprintf("Each simulation tests %d M-values.\n", length(M_values_to_test)))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", length(sac_levels) * n_sims * length(M_values_to_test) * n_reps))


##########
# 3. Define FULL Sites & Create Weight Matrix (w)
# (This is static geometry, so it lives OUTSIDE the loop)
##########

cat("Generating static full landscape geometry (1600 sites from 40000 cells)...\n")

# Assign each cell (1 to 40000) to a site (1 to 1600)
full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1
full_site_row <- (full_cell_row - 1) %/% site_dim + 1
full_site_col <- (full_cell_col - 1) %/% site_dim + 1

# This vector has length 40000, with values from 1 to 1600
full_site_id_for_cell <- (full_site_row - 1) * full_n_sites_x + full_site_col

# Create the full weight matrix w (1600 x 40000)
full_w <- matrix(0, full_M, full_n_cells)
for (i in 1:full_M) {
    full_w[i, full_site_id_for_cell == i] = 1
}

##########
# 4. Helper Function for Subsetting
##########

get_subset_landscape <- function(M_target, site_dim, full_cell_row, full_cell_col, full_cellCovs, full_site_id_for_cell, full_w, full_lambda_j) {
    
    # 1. Calculate dimensions of the target subset
    n_sites_per_dim_target <- sqrt(M_target)
    grid_dim_target <- n_sites_per_dim_target * site_dim
    
    # cat(sprintf("\n--- Subsetting for M = %d ---\n", M_target))
    
    # 2. Find which *full* cell indices are in the top-left quadrant
    cell_indices_to_keep <- which(full_cell_row <= grid_dim_target & full_cell_col <= grid_dim_target)
    
    # 3. Find which *full* site indices are in the top-left quadrant
    site_ids_in_subset <- unique(full_site_id_for_cell[cell_indices_to_keep])
    site_ids_to_keep <- sort(site_ids_in_subset) # Should be 1:M_target
    
    if (length(site_ids_to_keep) != M_target) {
        warning(sprintf("Expected %d sites, but found %d unique sites in subset.", M_target, length(site_ids_to_keep)))
    }
    
    # --- 4. Subset all the data ---
    
    # Subset cellCovs
    sub_cellCovs <- full_cellCovs[cell_indices_to_keep, , drop = FALSE]
    rownames(sub_cellCovs) <- NULL
    
    # Subset weight matrix 'w'
    sub_w <- full_w[site_ids_to_keep, cell_indices_to_keep]
    
    # Map site IDs
    global_site_ids_for_subset_cells <- full_site_id_for_cell[cell_indices_to_keep]
    local_site_id_map <- 1:M_target
    names(local_site_id_map) <- as.character(site_ids_to_keep)
    sub_site_id_for_cell <- as.numeric(local_site_id_map[as.character(global_site_ids_for_subset_cells)])
    
    # Subset cell coordinates
    sub_cell_row <- full_cell_row[cell_indices_to_keep]
    sub_cell_col <- full_cell_col[cell_indices_to_keep]
    
    # Subset lambda_j
    sub_lambda_j <- full_lambda_j[cell_indices_to_keep]
    
    return(list(
        M = M_target,
        n_cells = length(cell_indices_to_keep),
        grid_dim = grid_dim_target,
        n_sites_x = n_sites_per_dim_target,
        n_sites_y = n_sites_per_dim_target,
        cellCovs = sub_cellCovs,
        w = sub_w,
        site_id_for_cell = sub_site_id_for_cell,
        cell_row = sub_cell_row,
        cell_col = sub_cell_col,
        lambda_j = sub_lambda_j
    ))
}

##########
# 5. Initialize Loop & Storage
##########

all_results_df <- data.frame()

# Setup Output Directory
output_dir <- file.path("simulation_experiments", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n--- Starting Main Simulation Loop ---\n")

# --- Main loop over SAC Levels ---
for (sac_level in sac_levels) {
    
    current_sigma <- sac_sigmas[[sac_level]]
    all_plots_list <- list() # Store plots only for the current SAC level
    
    cat(sprintf("\n#######################################################\n"))
    cat(sprintf("### STARTING EXPERIMENT: SAC = %s (Sigma = %d) ###\n", sac_level, current_sigma))
    cat(sprintf("#######################################################\n"))

    # --- Loop over SIMULATIONS ---
    for (sim in 1:n_sims) {
        
        cat(sprintf("\n=== SAC: %s | Sim %d of %d ===\n", sac_level, sim, n_sims))
    
        ##########
        # 6. Create FULL Landscape (Cells & Covariates)
        #    WITH SPATIAL AUTOCORRELATION
        ##########
        
        # 1. Generate base random noise
        r_base <- terra::rast(nrows=full_grid_dim, ncols=full_grid_dim, 
                              xmin=0, xmax=full_grid_dim, ymin=0, ymax=full_grid_dim,
                              vals=rnorm(full_n_cells))
        
        # 2. Apply Gaussian smoothing if Sigma > 0 (Medium/High SAC)
        if (current_sigma > 0) {
            # Create Gaussian weights
            fw <- terra::focalMat(r_base, current_sigma, type = "Gauss")
            # Apply focal filter (this induces autocorrelation)
            r_smooth <- terra::focal(r_base, w = fw, fun = sum, na.rm = TRUE)
            
            # 3. Rescale to Mean=0, SD=1
            # (Important to keep effect sizes comparable across SAC levels)
            vals <- terra::values(r_smooth)
            vals_scaled <- as.vector(scale(vals))
            terra::values(r_smooth) <- vals_scaled
            
            r_final <- r_smooth
        } else {
            # Low SAC (Zero): Just ensure it's scaled (rnorm is already ~0,1)
            vals <- terra::values(r_base)
            terra::values(r_base) <- as.vector(scale(vals))
            r_final <- r_base
        }
        
        # Extract covariate vector (ensuring correct order)
        # Terra values are row-major (top-left to bottom-right), matching our grid indexing
        full_cellCovs <- data.frame(
            cell_cov1 = terra::values(r_final, mat=FALSE) 
        )
        # Handle edge-case NAs from focal (should be minimal/none with padding, but safer to fill)
        if(any(is.na(full_cellCovs$cell_cov1))) {
             full_cellCovs$cell_cov1[is.na(full_cellCovs$cell_cov1)] <- 0
        }
        
        # Create design matrix for state (X_design) for all cells
        full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
        
        # Calculate true latent abundance (lambda_j) for each cell
        full_log_lambda_j <- full_X_cell %*% true_betas
        full_lambda_j <- exp(full_log_lambda_j)
    
    
        # --- Main loop over M values ---
        for (M_i in M_values_to_test) {
            
            # --- 6.1. Get subset data for M = M_i ---
            subset_data <- get_subset_landscape(
                M_target = M_i,
                site_dim = site_dim,
                full_cell_row = full_cell_row,
                full_cell_col = full_cell_col,
                full_cellCovs = full_cellCovs,
                full_site_id_for_cell = full_site_id_for_cell,
                full_w = full_w,
                full_lambda_j = full_lambda_j
            )
            
            # --- 6.2. Unpack ---
            M <- subset_data$M
            w <- subset_data$w
            lambda_j <- subset_data$lambda_j
            cellCovs <- subset_data$cellCovs
            
            # cat(sprintf("  Running M=%d.\n", M))
            
            ##########
            # 7. Simulate True State (Occupancy Z)
            ##########
            
            lambda_tilde_i <- w %*% lambda_j
            psi_i <- 1 - exp(-lambda_tilde_i)
            Z_i <- rbinom(M, 1, psi_i)
            
            ##########
            # 8. Simulate Observation Data (y)
            ##########
            
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
                cellCovs = cellCovs,
                w = w
            )
            
            ##########
            # 10. Fit the occuN Model
            ##########
            
            best_fm <- NULL
            min_nll <- Inf
            n_params <- length(true_alphas) + length(true_betas)
        
            for (rep in 1:n_reps) {
                rand_starts <- runif(n_params, -5, 5) 
                fm_rep <- try(occuN(
                    formula = ~obs_cov1 ~ cell_cov1,
                    data = umf,
                    starts = rand_starts,
                    se = TRUE,
                    method = selected_optimizer
                ), silent = TRUE)
                
                if (inherits(fm_rep, "try-error")) next
                
                current_nll <- fm_rep@negLogLike
                if (current_nll < min_nll) {
                    min_nll <- current_nll
                    best_fm <- fm_rep
                }
            } 
            
            fm <- best_fm
            
            if (is.null(fm)) {
                    # cat(sprintf("    (All reps failed for M=%d)\n", M))
                    loop_results <- data.frame(
                        Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
                        True_Value = c(true_alphas, true_betas),
                        Estimated_Value = c(NA, NA, NA, NA),
                        M = M,
                        sim_id = sim,
                        SAC_Level = sac_level
                    )
                    all_results_df <- rbind(all_results_df, loop_results)
                    
                    if(sim == 1) all_plots_list <- c(all_plots_list, list(ggplot(), ggplot(), ggplot()))
                    next 
            }
        
            ##########
            # 11. Store Results
            ##########
            
            loop_results <- data.frame(
                Parameter = c("alpha (det_int)", "alpha (det_cov1)", "beta (state_int)", "beta (state_cov1)"),
                True_Value = c(true_alphas, true_betas),
                Estimated_Value = c(coef(fm, 'det'), coef(fm, 'state')),
                M = M, 
                sim_id = sim,
                SAC_Level = sac_level # NEW: Store SAC Level
            )
            
            all_results_df <- rbind(all_results_df, loop_results)
            
            ##########
            # 12. Generate Plots (ONLY FOR SIM 1 of CURRENT SAC)
            ##########
            
            if (sim == 1) {
                # Create data frame for cell covariates
                cell_df <- data.frame(
                    x = subset_data$cell_col,
                    y = subset_data$cell_row,
                    covariate = cellCovs$cell_cov1
                )
                
                # Map site-level results
                cell_df$site_latent_abundance <- as.numeric(lambda_tilde_i[subset_data$site_id_for_cell])
                cell_df$site_occupancy_prob <- as.numeric(psi_i[subset_data$site_id_for_cell])
                cell_df$site_true_occupancy <- as.factor(Z_i[subset_data$site_id_for_cell])
                
                # Site boxes
                xmin_vals <- seq(from = 0.5, to = subset_data$grid_dim - site_dim + 0.5, by = site_dim)
                xmax_vals <- seq(from = site_dim + 0.5, to = subset_data$grid_dim + 0.5, by = site_dim)
                site_boxes <- expand.grid(xmin = xmin_vals, ymin = xmin_vals)
                site_boxes$xmax <- site_boxes$xmin + site_dim
                site_boxes$ymax <- site_boxes$ymin + site_dim
                
                # Plot 1: Covariate
                p_covariate <- ggplot(cell_df, aes(x = x, y = y, fill = covariate)) +
                    geom_raster() +
                    scale_fill_viridis_c() +
                    coord_fixed(expand = FALSE) +
                    geom_rect(data = site_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                        color = "red", fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
                    labs(title = sprintf("Sites & Covariate (M=%d, SAC=%s)", M, sac_level), fill = "Covariate") +
                    theme_minimal()
                
                # Plot 2: Abundance
                p_abundance <- ggplot(cell_df, aes(x = x, y = y, fill = site_latent_abundance)) +
                    geom_raster() +
                    scale_fill_viridis_c(option = "magma") +
                    coord_fixed(expand = FALSE) +
                    geom_rect(data = site_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                        color = "red", fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
                    labs(title = sprintf("Site Abundance (M=%d)", M), fill = "Latent Abund.") +
                    theme_minimal()
                
                # Plot 3: Occupancy
                p_occupancy <- ggplot(cell_df, aes(x = x, y = y, fill = site_true_occupancy)) +
                    geom_raster() +
                    scale_fill_manual(values = c("0" = "navyblue", "1" = "yellow")) +
                    coord_fixed(expand = FALSE) +
                    geom_rect(data = site_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                        color = "red", fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
                    labs(title = sprintf("Site Occupancy (M=%d)", M), fill = "Occupied") +
                    theme_minimal()
                
                all_plots_list <- c(all_plots_list, list(p_covariate, p_abundance, p_occupancy))
            } 
            
        } # --- End M Loop ---
        
    } # --- End Sim Loop ---
    
    ##########
    # 13. Save SAC-Specific Landscape Plot
    ##########
    
    cat(sprintf("\nSaving landscape plots for SAC=%s...\n", sac_level))
    combined_plot <- patchwork::wrap_plots(all_plots_list, nrow = length(M_values_to_test), ncol = 3)
    
    ggsave(file.path(output_dir, sprintf("plot_SAC=%s.png", sac_level)), 
                 plot = combined_plot, dpi = 300, width = 18, height = 26) 
    
} # --- End SAC Level Loop ---


##########
# 14. Save Aggregate Results & Error Boxplots
##########

cat("\n--- Simulation Study Complete ---\n")

# Save the full results data frame
write.csv(all_results_df, file.path(output_dir, "params.csv"), row.names = FALSE)
cat(sprintf("All parameters saved to %s/params.csv\n", output_dir))

cat("Generating combined error boxplots (colored by SAC level)...\n")

# Calculate error
all_results_df$Error <- all_results_df$True_Value - all_results_df$Estimated_Value
all_results_df$M_factor <- as.factor(all_results_df$M)
# Ensure SAC factor order
all_results_df$SAC_Level <- factor(all_results_df$SAC_Level, levels = c("Low", "Medium", "High"))

# Helper for boxplots
create_error_plot <- function(param_name, title) {
  ggplot(all_results_df[all_results_df$Parameter == param_name, ], 
         aes(x = M_factor, y = Error, fill = SAC_Level)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    scale_fill_manual(values = c("Low" = "yellow", "Medium" = "orange", "High" = "red")) +
    labs(title = title, x = "M (Sites)", y = "Error (True - Est.)") +
    theme_bw()
}

p_err_beta_int <- create_error_plot("beta (state_int)", "State Intercept")
p_err_beta_cov <- create_error_plot("beta (state_cov1)", "State Slope")
p_err_alpha_int <- create_error_plot("alpha (det_int)", "Observation Intercept")
p_err_alpha_cov <- create_error_plot("alpha (det_cov1)", "Observation Slope")

# Combine the 4 error plots
combined_error_plot <- (p_err_beta_int | p_err_beta_cov) / (p_err_alpha_int | p_err_alpha_cov) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save the combined error plot
ggsave(file.path(output_dir, "error_boxplots.png"), 
             plot = combined_error_plot, dpi = 300, width = 12, height = 12)

cat(sprintf("Error boxplots saved to %s/error_boxplots.png\n", output_dir))
cat("--- Script Finished ---\n")
