####################
# Helper functions

# November 9, 2025
####################

library(terra) # geospatial operations


#######
# extract environmental features
# at checklist locations
#######
extract_state_covs <- function(df, cov_tif, x = "longitude", y = "latitude", crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") {
  # Convert the dataframe to a SpatVector
  df_pts <- vect(df, geom = c(x, y), crs = crs)
  
  # Extract environmental features
  env_vars_df <- data.frame(
    checklist_id = df$checklist_id,
    terra::extract(cov_tif, df_pts)
  )
  return(env_vars_df)
}

#######
# normalize dataset
#######
norm_ds <- function(df, obs_covs, state_covs, norm_list = list()){

  if(length(norm_list) == 0){
    for(name in c(obs_covs, state_covs)){
      # calc mean/var for each cov, if training
      ma <- max(df[[name]])
      mi <- min(df[[name]])
      norm_list[[name]] <- c(ma, mi)
    }
  }
  
  # xi - min(x)/(max(x) - min(x))
  for(cov in names(norm_list)){
    df[[cov]] <- (df[[cov]] - norm_list[[cov]][[2]])/(norm_list[[cov]][[1]] - norm_list[[cov]][[2]])
  }
  
  return(list(df=df, n_l=norm_list))
  
}

#######
# spatial subsampling as defined by: 
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/ddi.13271
#######
spatial_subsample <- function(df, cell_names){
  valid_df <- data.frame()
  i <- 0
  for(freq in table(df$cell)){ 
    i <- i + 1
    if(freq > 1){
      checklists <- df[df$cell == cell_names[i],]
      sample <- checklists[sample(nrow(checklists), 1), ]
    } else {
      sample <- df[df$cell == cell_names[i],]
    }
    valid_df <- rbind(valid_df, sample)
  }  
  return(valid_df)
}




#########
# Rounding Lat/Long
#########
round_lat_long <- function(df, rounding_degree){
  df$rounded_lat <- round(df$latitude, digits = rounding_degree)
  df$rounded_long <- round(df$longitude, digits = rounding_degree)
  df$rounded_locality_id <- paste(as.character(df$rounded_long), as.character(df$rounded_lat), sep = "_")
  
  # Remove the temporary columns by setting them to NULL
  df$rounded_lat <- NULL
  df$rounded_long <- NULL
  
  return(df)
}