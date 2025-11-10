
library(dplyr)

source(file.path("R","utils.R"))



prepare_train_data <- function (cov_tif, state_covs, obs_covs, placeholder_spec_name = "AMCR"){

  train_filename <- paste0(placeholder_spec_name, "_zf_filtered_region_2017.csv")
  train_df_og <- read.delim(file.path("checklist_data","species", placeholder_spec_name, train_filename), sep = ",", header = T)
  
  train_df_og <- train_df_og[!is.na(train_df_og$duration_minutes),]

  train_df_og <- train_df_og[train_df_og$observation_date >= "2017-05-15" & train_df_og$observation_date <= "2017-07-09",]
   
  train_df <- train_df_og


  train_env_df <- extract_state_covs(train_df, cov_tif)
  train_df <- inner_join(train_df, train_env_df, by = "checklist_id")
  norm_res <- norm_ds(train_df, obs_covs, state_covs)

  train_df_unnorm <- train_df[,c(obs_covs, state_covs)]

  train_df <- norm_res$df
  norm_list <- norm_res$n_l
  
  for(name in c(obs_covs, state_covs)){
    train_df[,paste("unnorm_",name)] = train_df_unnorm[,name]
    
  }

  train_df$species_observed <- -1
  train_df$occupied_prob <- -1
  train_df$det_prob <- -1
  
   
  train_df$formatted_date <- train_df$observation_date

  return (list(train_df = train_df, norm_list = norm_list))
}



simulate_train_data <-  function (sites_df, occ_par_list, det_par_list){
  
  sites_df_u <- subset(sites_df, !duplicated(site))
  sites_list <- sites_df_u$site

  sites_df_u$occupied_prob <- calculate_weighted_sum(occ_par_list, sites_df_u)

  sites_df_u$occupied_prob <- expit(sites_df_u$occupied_prob)
  sites_df_u$occupied <- rbinom(nrow(sites_df_u), 1, sites_df_u$occupied_prob)

  j<-0

  res_df <- NA
  for(c_site in sites_list){
    
    j = j+1

    checklists_at_site <- sites_df[sites_df$site == c_site,]
  
    checklists_at_site$occ_prob <- sites_df_u[sites_df_u$site == c_site,]$occ_prob
    checklists_at_site$occupied <- sites_df_u[sites_df_u$site == c_site,]$occupied
    
    if(j==1){
      res_df = checklists_at_site
    } else {
      res_df = rbind(res_df, checklists_at_site)
    }
  
  }


  res_df$det_prob <- calculate_weighted_sum(det_par_list, res_df)
    
  res_df$det_prob <- expit(res_df$det_prob)

  res_df$detection <- rbinom(nrow(res_df), 1, res_df$det_prob)

  res_df$species_observed <- res_df$occupied * res_df$detection


  
  cov_names <- c(occ_covs, det_covs)

  res_df <- res_df[,!names(res_df) %in% cov_names]


  for(name in c(det_covs, occ_covs)){
    res_df[,name] = res_df[,paste("unnorm_",name)]
    
  }

  return (res_df)
}

simulate_test_data <- function (norm.list, occ_par_list, det_par_list, placeholder_spec_name = "AMCR"){
  
  f.name_test <- paste0(placeholder_spec_name, "/", placeholder_spec_name, "_zf_filtered_region_2018.csv")

  test.df.og <- read.delim(paste0("checklist_data/species/", f.name_test), sep = ",", header = T)
  
  test.df.og <- test.df.og[!is.na(test.df.og$duration_minutes),]

  test.df.og <- test.df.og[test.df.og$observation_date >= "2018-05-15" & test.df.og$observation_date <= "2018-07-09",]

  test.df <- test.df.og

  
  test_env.df <- extractEnvFeat(test.df, OR.init, obs_covs = det_covs)
  test.df <- inner_join(test.df, test_env.df, by = "checklist_id")
  norm_res <- norm_ds(test.df, det_covs, occ_covs, test=TRUE, norm.list = norm.list)
  
  test.df_unnorm <- test.df[,c(det_covs, occ_covs)]

  test.df <- norm_res$df

  for(name in c(det_covs, occ_covs)){
    test.df[,paste("unnorm_",name)] = test.df_unnorm[,name]
    
  }
  
  test.df$species_observed <- -1
  test.df$occupied_prob <- -1
  test.df$det_prob <- -1
  

  test.df$occupied_prob <- calculate_weighted_sum(occ_par_list, test.df)

  test.df$occupied_prob <- expit(test.df$occupied_prob)
  test.df$occupied <- rbinom(nrow(test.df), 1, test.df$occupied_prob)


  test.df$det_prob <- calculate_weighted_sum(det_par_list, test.df)
    
  test.df$det_prob <- expit(test.df$det_prob)

  test.df$detection <- rbinom(nrow(test.df), 1, test.df$det_prob)

  test.df$species_observed <- test.df$occupied * test.df$detection

  
  cov_names <- c(occ_covs, det_covs)

  test.df <- test.df[,!names(test.df) %in% cov_names]

  for(name in c(det_covs, occ_covs)){
    test.df[,name] = test.df[,paste("unnorm_",name)]
    
  }

  return (test.df)
}
