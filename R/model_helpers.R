


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


