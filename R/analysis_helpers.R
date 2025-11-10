

calcClusteringStats <- function(pred_df, og_df){
    
  
    pred_df <- pred_df[order(pred_df$checklist_id),]
    og_df <- og_df[order(og_df$checklist_id),]
    
    pred_sites <- as.factor(pred_df$site)
    og_sites <- as.factor(og_df$site)

    ari <- adjustedRandIndex(og_sites, pred_sites)
    ami <- AMI(og_sites, pred_sites)
    nid <- NID(og_sites, pred_sites)
    
 
    return(list(ari=ari, ami=ami , nid=nid))

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

