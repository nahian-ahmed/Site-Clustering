
source(file.path("R","utils.R"))



prep_train_data <- function (cov_tif, placeholder_spec_name = "AMCR"){

    f.name_train <- paste0(placeholder_spec_name, "/", placeholder_spec_name, "_zf_filtered_region_2017.csv")

    train.df.og <- read.delim(paste0("checklist_data/species/", f.name_train), sep = ",", header = T)
  
    train.df.og <- train.df.og[!is.na(train.df.og$duration_minutes),]

    train.df.og <- train.df.og[train.df.og$observation_date >= "2017-05-15" & train.df.og$observation_date <= "2017-07-09",]
   
    train.df <- train.df.og


    train_env.df <- extract_state_covs(train.df, cov_tif)
    train.df <- inner_join(train.df, train_env.df, by = "checklist_id")
    norm_res <- norm_ds(train.df, det_covs, occ_covs)

    train.df_unnorm <- train.df[,c(det_covs, occ_covs)]

    train.df <- norm_res$df
    norm.list <- norm_res$n_l
    
    for(name in c(det_covs, occ_covs)){
        train.df[,paste("unnorm_",name)] = train.df_unnorm[,name]
        
    }

    train.df$species_observed <- -1
    train.df$occupied_prob <- -1
    train.df$det_prob <- -1
    
   
    train.df$formatted_date <- train.df$observation_date

    return (list(train.df = train.df, norm.list = norm.list))
}