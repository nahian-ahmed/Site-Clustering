####################
# Helper functions

# April 18, 2025
####################

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(sf) # spatial geometry operations
library(auk) # ebird data processing
library(unmarked) # occupancy modeling
library(rje) # expit function
library(terra) # geospatial operations
library(mclust) # ARI calculation
library(aricode) # AMI and NID calculation



#########
# Rounding Lat/Long
#########
roundLatLong <- function(df, rounding_degree){
    df$rounded_lat <- round(df$latitude, digits = rounding_degree)
    df$rounded_long <- round(df$longitude, digits = rounding_degree)
    df$rounded_locality_id <- paste(as.character(df$rounded_long), as.character(df$rounded_lat), sep = "_")
    
    # Remove the temporary columns by setting them to NULL
    df$rounded_lat <- NULL
    df$rounded_long <- NULL
    
    return(df)
}


##########
# site closure for Occ Model:
#     1. constant site covariates
#     2. no false positives (detected only if occupied)
##########
enforceClosure <- function(sites_df, occ_cov_list, sites_list){
    j<-1
    closed_df <- NA

    for(eBird_site in sites_list){
        
    
        checklists_at_site <- sites_df[sites_df$site == eBird_site,]
        
        for(occCov_i in occ_cov_list){
            checklists_at_site[occCov_i] <- mean(checklists_at_site[[occCov_i]])
        
        }

        
        if(j==1){
            closed_df = checklists_at_site
        } else {
            closed_df = rbind(closed_df, checklists_at_site)
        }
        j = j+1
    
    }
    return(closed_df)
}

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



#######
# extract environmental features
# at checklist locations
#######
extractEnvFeat <- function(df, OR.tif, obs_covs) {
    # Convert the dataframe to a SpatVector
    df.pts <- vect(df, geom = c("longitude", "latitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    # Extract environmental features
    env_vars.df <- data.frame(
        checklist_id = df$checklist_id,
        terra::extract(OR.tif, df.pts)
    )
    return(env_vars.df)
}


#######
# normalize dataset
#######
norm_ds <- function(df, det_covs, occ_covs, test=FALSE, norm.list=list()){

    if(length(norm.list) == 0){
        for(name in c(det_covs, occ_covs)){
            # calc mean/var for each cov, if training
            ma <- max(df[[name]])
            mi <- min(df[[name]])
            norm.list[[name]] <- c(ma, mi)
        }
    }
    
    # xi - min(x)/(max(x) - min(x))
    for(cov in names(norm.list)){
        df[[cov]] <- (df[[cov]] - norm.list[[cov]][[2]])/(norm.list[[cov]][[1]] - norm.list[[cov]][[2]])
    }
    
    return(list(df=df, n_l=norm.list))
    
}

#######
# spatial subsampling as defined by: 
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/ddi.13271
#######
spatial.subsample <- function(df, cell.names){
    valid.df <- data.frame()
    i <- 0
    for(freq in table(df$cell)){ 
        i <- i + 1
        if(freq > 1){
            chklsts <- df[df$cell == cell.names[i],]
            samp <- chklsts[sample(seq(1:nrow(chklsts)), 1),]
        } else {
            samp <- df[df$cell == cell.names[i],]
        }
        valid.df <- rbind(valid.df, samp)
    }    
    return(valid.df)
}



#######
# calculates the occupancy model from a given dataset
# containing checklists and sites
#######
# 1. enforces closure
# 2. formats it w/r/t eBird data
# 3. runs through occupancy model
#######
calcOccModel <- function(df, occ_covs, det_covs, skip_closure=FALSE){

    sites_occ <- subset(df, !duplicated(site))$site

    
    closed_df <- df
    if(!skip_closure){
        closed_df <- enforceClosure(df, occ_covs, sites_occ)
    } 
    
    
    umf_AUK <- auk::format_unmarked_occu(
        closed_df,
        site_id = "site",
        response = "species_observed",
        site_covs = occ_covs,
        obs_covs = det_covs
    )
    
    det_cov_str <- paste("", paste(det_covs, collapse="+"), sep=" ~ ")
    occ_cov_str <- paste("", paste(occ_covs, collapse="+"), sep=" ~ ")
    
    species_formula <- paste(det_cov_str, occ_cov_str, sep = " ")
    species_formula <- as.formula(species_formula)
    
    occ_um <- unmarked::formatWide(umf_AUK, type = "unmarkedFrameOccu")

    og_syn_gen_form <- unmarked::occu(formula = species_formula, data = occ_um)
  
    return(og_syn_gen_form)

}
#######





########
predict_sdm_map <- function(occ_pars, region, intercept = TRUE){
    
    valid_boundary <- terra::vect("occupancy_feature_raster/boundary/boundary.shp")
    crs(valid_boundary) <- crs(region)
    region <- terra::crop(region, valid_boundary, mask = TRUE)
    
    weighted_sum <- 0
    par_idx <- 1
    n_pars <- length(occ_pars)

    if(intercept){
        weighted_sum <- weighted_sum + occ_pars[[1]]
        par_idx <- par_idx + 1
    }
    while(par_idx <= n_pars){
        weighted_sum <- weighted_sum + (occ_pars[[par_idx]] * subset(region, names(occ_pars)[[par_idx]]))
        par_idx <- par_idx + 1
    }

    region$occ_prob <- expit(weighted_sum)


    return(region)
}
########

########
get_occu_map_diff <- function(occu_map_gt, occu_map) {

    occu_m_gt <- occu_map_gt[["occu"]]
    occu_m <- occu_map[["occu"]]

    difference <- occu_m - occu_m_gt
    percentage_difference <- (abs(occu_m - occu_m_gt) / occu_m_gt) * 100

    # Create a raster stack with two bands
    occu_map_diff <- rast(nrows = nrow(occu_m_gt), ncols = ncol(occu_m_gt), 
                          ext = ext(occu_m_gt), crs = crs(occu_m_gt), nlyrs = 2)
    
    values(occu_map_diff)[,1] <- values(difference)
    values(occu_map_diff)[,2] <- values(percentage_difference)

    names(occu_map_diff) <- c("difference", "percentage_difference")

    return(occu_map_diff)
}
########


########
calculate_weighted_sum <- function(pars, covs_df, intercept = TRUE){

    weighted_sum <- 0
    par_idx <- 1
    n_pars <- length(pars)

    if(intercept){
        weighted_sum <- weighted_sum + pars[[1]]
        par_idx <- par_idx + 1
    }
    while(par_idx <= n_pars){
        weighted_sum <- weighted_sum + (pars[[par_idx]] * covs_df[[names(pars)[[par_idx]]]])
        par_idx <- par_idx + 1
    }
    return(weighted_sum)
}


########
get_parameters <- function(df, i, occ_covs, det_covs, occ_intercept = TRUE, det_intercept = TRUE){

    occ_par_list <- list()
    if (occ_intercept){
        occ_par_list[["occ_intercept"]] <- df[i, "occ_intercept"] 
    }
    for (occ_cov in occ_covs){
        occ_par_list[[occ_cov]] <- df[i, occ_cov]

    }

    det_par_list <- list()
    if (det_intercept){
        det_par_list[["det_intercept"]] <- df[i, "det_intercept"] 
    }
    for (det_cov in det_covs){
        det_par_list[[det_cov]] <- df[i, det_cov]

    }

 
    return (list(occ_par_list = occ_par_list, det_par_list = det_par_list))
}
########













