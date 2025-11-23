library(aricode)
library(PRROC)
library(dplyr)
library(sf)



calculate_clustering_stats <- function(ref_df, comp_df) {
  
  # Ensure we have the required libraries
  if (!requireNamespace("aricode", quietly = TRUE)) {
    warning("Package 'aricode' is required for AMI and NID. Returning NAs.")
    return(list(ARI = NA, AMI = NA, NID = NA))
  }
  
  # 1. Align the dataframes by checklist_id
  # We use inner_join to ensure we are comparing the exact same checklists
  # Renaming columns to avoid collision
  aligned_data <- dplyr::inner_join(
    ref_df %>% dplyr::select(checklist_id, ref_label = site),
    comp_df %>% dplyr::select(checklist_id, comp_label = site),
    by = "checklist_id"
  )
  
  if (nrow(aligned_data) == 0) {
    warning("No matching checklists found between the two clusterings.")
    return(list(ARI = NA, AMI = NA, NID = NA))
  }
  
  # 2. Calculate Metrics
  c1 <- aligned_data$ref_label
  c2 <- aligned_data$comp_label
  
  # Adjusted Rand Index
  ari_val <- aricode::ARI(c1, c2)
  
  # Adjusted Mutual Information
  ami_val <- aricode::AMI(c1, c2)
  
  # Normalized Information Distance (1 - NMI)
  # aricode::NID usually calculates 1 - NMI. 
  nid_val <- aricode::NID(c1, c2)
  
  return(list(
    ARI = ari_val,
    AMI = ami_val,
    NID = nid_val
  ))
}




#' Calculate Classification Metrics (AUC/AUPRC)
#'
#' Calculates AUC-ROC and AUC-PR based on predicted probabilities vs truth.
#' Returns NA if calculation fails or returns empty values.
calculate_classification_metrics <- function(pred_prob, true_labels) {
  
  # 1. Basic length validation
  if (length(pred_prob) != length(true_labels)) {
    warning("Length of predictions and truth do not match")
    return(list(auc = NA, auprc = NA))
  }
  
  # 2. Handle NAs in predictions (treat as 0 or skip)
  # Removing NAs ensures PRROC doesn't crash
  valid_idx <- !is.na(pred_prob) & !is.na(true_labels)
  pred_prob <- pred_prob[valid_idx]
  true_labels <- true_labels[valid_idx]

  # 3. Handle edge case: Only 1 class present (e.g., all 0s)
  if (length(unique(true_labels)) < 2) {
    return(list(auc = NA, auprc = NA))
  }
  
  # 4. Attempt PRROC calculation
  pr_metrics <- try({
    PRROC::pr.curve(
      scores.class0 = pred_prob[true_labels == 1], # Positives
      scores.class1 = pred_prob[true_labels == 0], # Negatives
      curve = FALSE
    )
  }, silent = TRUE)
  
  # 5. Extract values safely
  # If pr_metrics is error OR if the specific slot is NULL, return NA
  
  auc_val <- NA
  if (!inherits(pr_metrics, "try-error") && !is.null(pr_metrics$auc.roc)) {
    auc_val <- pr_metrics$auc.roc
  }

  auprc_val <- NA
  if (!inherits(pr_metrics, "try-error") && !is.null(pr_metrics$auc.integral)) {
    auprc_val <- pr_metrics$auc.integral
  }
  
  return(list(
    auc = auc_val,
    auprc = auprc_val
  ))
}


calcDescriptiveClusteringStatsWithReference <- function(df, cluster_col, feature_cols, normalize = FALSE) {
    
    if (normalize) {
        df[feature_cols] <- lapply(df[feature_cols], function(x) {
            (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        })
    }
    
    cluster_sizes <- df %>%
        group_by(!!as.name(cluster_col)) %>%
        summarise(freq = n(), .groups = "drop")



 
    descr_stats <- list(
        n_points = nrow(df),
        n_clusters = nrow(cluster_sizes),
        min_size = min(cluster_sizes$freq),
        max_size = max(cluster_sizes$freq),
        mean_size = round(mean(cluster_sizes$freq), 4),
        sd_size = round(sd(cluster_sizes$freq), 4),
        perc_svs = round((sum(cluster_sizes$freq == 1) / nrow(cluster_sizes)) * 100, 4)
    )
    
    return(descr_stats)
}

calcDescriptiveClusteringStats <- function(clustered_df) {
    
        # Calculate clustering and species specific stats
    clust_freq_df <- clustered_df %>% group_by(site) %>% dplyr::summarise(freq=n())

    descr_stats <- list(
        n_points = nrow(clustered_df),
        n_clusters = nrow(clust_freq_df),
        min_size = min(clust_freq_df$freq),
        max_size = max(clust_freq_df$freq),
        mean_size = round(mean(clust_freq_df$freq), 4),
        sd_size = round(sd(clust_freq_df$freq), 4),
        perc_svs = round((nrow(clust_freq_df[clust_freq_df$freq==1,])/nrow(clust_freq_df)) *100, 4)
     )
    return(descr_stats)
}



summarize_clusterings <- function(all_clusterings, all_site_geometries, units = "m") {
  
  all_summaries_list <- list()
  
  # Set unit divisors
  area_divisor <- 1.0
  area_unit_label <- "m2"
  
  if (units == "km") {
    area_divisor <- 1000000.0
    area_unit_label <- "km2"
  }
  
  for (method_name in names(all_clusterings)) {
    
    cat(paste("--- Summarizing clustering for:", method_name, "---\n"))
    
    # 1. Get the correct points dataframe
    cluster_data <- all_clusterings[[method_name]]
    if (is.list(cluster_data) && "result_df" %in% names(cluster_data)) {
      cluster_data <- cluster_data$result_df
    }
    
    # 2. Get the correct geometries dataframe
    geom_data <- all_site_geometries[[method_name]]
    
    # Skip if either is missing (e.g., a method failed)
    if (is.null(cluster_data) || is.null(geom_data)) {
      cat(paste("    - WARNING: Missing cluster or geometry data for", method_name, ". Skipping.\n"))
      next
    }
    
    # --- 3. Calculate Point Stats (Metrics 1-5) ---
    # (n_points per site)
    point_stats <- cluster_data %>%
      group_by(site) %>%
      dplyr::summarise(n_points = n())
      
 
    # --- 5. Calculate Area Stats (Metrics 10-13) ---
    # (area of the site geometry)
    # Geometries are in Albers (meters), so st_area() returns m^2
    geom_data$area_m2 <- as.numeric(sf::st_area(geom_data))
    
    # --- 6. Join all stats by site ---
    # We use geom_data as the base, as it represents the definitive list of sites
    site_summary <- geom_data %>%
      dplyr::left_join(point_stats, by = "site") 
      
    # Handle sites that might have had 0 points (n_points=NA)
    # site_summary$n_points[is.na(site_summary$n_points)] <- 0

      
    # --- 7. Apply Unit Conversions ---
    site_summary$area <- site_summary$area_m2 / area_divisor
    
    # --- 8. Aggregate final metrics for the method ---
    method_summary <- site_summary %>%
      sf::st_drop_geometry() %>% # Drop geometry before summarizing
      dplyr::summarise(
        n_sites = n(),
        min_points = min(n_points),
        max_points = max(n_points),
        mean_points = mean(n_points),
        
        min_area = min(area),
        max_area = max(area),
        mean_area = mean(area)
      )
    
    method_summary$method <- method_name
    
    # Add to our list
    all_summaries_list[[method_name]] <- method_summary
  }
  
  # Combine all method summaries into one dataframe
  final_df <- dplyr::bind_rows(all_summaries_list)
  
  # Rename columns to include units
  colnames(final_df) <- gsub("_area", paste0("_area_", area_unit_label), colnames(final_df))
  
  # Reorder columns to be logical
  final_df <- final_df %>%
    dplyr::select(
      method,
      n_sites,
      min_points,
      max_points,
      mean_points,
      dplyr::starts_with("min_area"),
      dplyr::starts_with("max_area"),
      dplyr::starts_with("mean_area"),
    )
  
  return(final_df)
}




#' Calculate Descriptive Statistics for Simulated Datasets
#'
#' This function takes simulated training and test dataframes and calculates
#' key descriptive statistics for both, returning them as a single-row dataframe.
#'
#' @param train_data The dataframe of simulated training checklists,
#'   as generated by `simulate_train_data`. Must contain 'site',
#'   'occupied', 'detection', 'species_observed', and 'N'.
#' @param test_data The dataframe of simulated test checklists,
#'   as generated by `simulate_test_data`. Must contain
#'   'occupied', 'detection', 'species_observed', and 'N'.
#'
#' @return A 1-row data.frame containing the 8 calculated statistics.
summarize_datasets <- function(train_data, test_data) {
  
  # --- 1. Calculate Training Statistics ---
  
  # Get unique site-level data from the training set.
  # 'occupied' (Z_i) and 'N' (N_i) are constant for all
  # checklists at a given site.
  train_sites <- dplyr::distinct(train_data, site, .keep_all = TRUE)
  
  # 1. Mean site-level abundance (N_i)
  train_site_mean_abundance <- mean(train_sites$N, na.rm = TRUE)
  
  # 2. Site-level occupancy rate (Z_i)
  train_site_occupancy_rate <- mean(train_sites$occupied, na.rm = TRUE)
  
  # 3. Point-level detection rate (d_it)
  train_point_detection_rate <- mean(train_data$detection, na.rm = TRUE)
  
  # 4. Point-level prevalence (y_it)
  train_point_prevalence <- mean(train_data$species_observed, na.rm = TRUE)
  
  
  # --- 2. Calculate Test Statistics ---
  # For test data, each point is effectively its own "site".
  
  # 5. Mean point-level abundance (N_j)
  test_point_mean_abundance <- mean(test_data$N, na.rm = TRUE)
  
  # 6. Point-level occupancy rate (Z_j)
  test_point_occupancy_rate <- mean(test_data$occupied, na.rm = TRUE)
  
  # 7. Point-level detection rate (d_j)
  test_point_detection_rate <- mean(test_data$detection, na.rm = TRUE)
  
  # 8. Point-level prevalence (y_j)
  test_point_prevalence <- mean(test_data$species_observed, na.rm = TRUE)
  
  
  # --- 3. Combine and Return ---
  
  # Return as a single-row data.frame
  stats_row <- data.frame(
    train_site_mean_abundance = train_site_mean_abundance,
    train_site_occupancy_rate = train_site_occupancy_rate,
    train_point_detection_rate = train_point_detection_rate,
    train_point_prevalence = train_point_prevalence,
    test_point_mean_abundance = test_point_mean_abundance,
    test_point_occupancy_rate = test_point_occupancy_rate,
    test_point_detection_rate = test_point_detection_rate,
    test_point_prevalence = test_point_prevalence
  )
  
  return(stats_row)
}