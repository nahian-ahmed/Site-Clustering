#############################################################################
# File to run spatial clustering algorithms on semi-simulated species data

# February 20, 2025
############################################################################



library(rje) # expit function
library(PRROC) # AUC and AUPRC calculation

species_df <- read.delim("simulation_species.csv", sep = ",", header = T)


train.df.og <- read.delim(paste0("results/simulation_experiments/raw/data/", spec_name, "_train.csv"), sep = ",", header = T)


method_names <- as.list(method_names)


methods <- genExp(method_names)


test_split_path = "results/simulation_experiments/raw/test_data_splits/"


model_list <- list()

for(exp_run in seq(1:runs)){

    train.df <- train.df.og
    test.df <- read.delim(paste0(test_split_path, "/", spec_name, "/run=", exp_run,".csv"), sep = ",", header = T)
    
    cat(paste0("\n\n\n", spec_name, " running experiment #", exp_run, "\n\n"))
    
  
    

    
    norm_res <- norm_ds(train.df, det_covs, occ_covs)
    train.df <- norm_res$df
    norm.list <- norm_res$n_l
    
    
    MDD <- 250
    MAX_OBS <- 100000
    MIN_OBS <- 1

    norm_res <- norm_ds(test.df, det_covs, occ_covs, test=TRUE, norm.list = norm.list)
    test.df <- norm_res$df

    train.df.gt <- train.df
    train.df <- subset(train.df, select = -c(site))
    
    if(exp_run == 1){
        groupedSite <- getClusterings(methods, train.df, occ_covs, det_covs, train.df.gt)
    }

    
    for(method_name in names(groupedSite)){
        

        set.seed(1)



        if(exp_run == 1){
            
            
        
            if (startsWith(method_name, "BayesOpt")){

                best_par_df <- data.frame(parameter = names(groupedSite[[method_name]]$Best_Pars), value = as.numeric(groupedSite[[method_name]]$Best_Pars))
                
                write.csv(best_par_df, paste0("results/simulation_experiments/raw/clustering_parameters/", spec_name,"_", method_name,".csv"), row.names=FALSE)

                groupedSite[[method_name]] <- groupedSite[[method_name]]$result_df 

            }
            if (!(method_name %in% c("1to10", "2to10", "2to10-sameObs", "1-UL"))){
                
                cl_stats <- calcClusteringStats(groupedSite[[method_name]], train.df.og)
            }
            else {
                cl_stats <- list(ari=NA, ami=NA , nid=NA)
            }

            cl_stats <- c(cl_stats, calcDescriptiveClusteringStatsWithReference(groupedSite[[method_name]], "site", occ_covs, normalize = FALSE))
            write.csv(data.frame(cl_stats), paste0("results/simulation_experiments/raw/clusterings/", spec_name,"_", method_name,"_clust_stats.csv"), row.names=FALSE)
            write.csv(groupedSite[[method_name]][,c("checklist_id", "longitude", "latitude", "site")], paste0("results/simulation_experiments/raw/clusterings/", spec_name,"_", method_name,".csv"), row.names=FALSE)

            model_list[[method_name]] <- list()
            
            test.formula <- calcOccModel(groupedSite[[method_name]], occ_covs, det_covs)

            occ_par_list <- test.formula@estimates@estimates$state@estimates 
            det_par_list <- test.formula@estimates@estimates$det@estimates

            model_list[[method_name]][["occu_parameters"]] <- occ_par_list
            model_list[[method_name]][["det_parameters"]] <- det_par_list 
            
            occ_par_mape <- mean(abs((unlist(occ_par_list[1:(1+length(occ_covs))]) - as.numeric(species_df[species_df$species==spec_name, c("occ_intercept", occ_covs)]))/as.numeric(species_df[species_df$species==spec_name, c("occ_intercept", occ_covs)]))) * 100
            det_par_mape <- mean(abs((unlist(det_par_list[1:(1+length(det_covs))]) - as.numeric(species_df[species_df$species==spec_name, c("det_intercept", det_covs)]))/as.numeric(species_df[species_df$species==spec_name, c("det_intercept", det_covs)]))) * 100

            model_pars_df <- data.frame(
                species = spec_name,
                method = method_name
            )
    
            model_pars_df$occ_intercept <- occ_par_list[[1]]
            par_idx <- 2
            for (occ_cov in occ_covs){
                model_pars_df[, occ_cov] <- occ_par_list[[par_idx]] 
                par_idx <- par_idx + 1
            }

            model_pars_df$det_intercept <- det_par_list[[1]]
            par_idx <- 2
            for (det_cov in det_covs){
                model_pars_df[, det_cov] <- det_par_list[[par_idx]] 
                par_idx <- par_idx + 1
            }
             
            model_pars_df$occ_par_mape <- occ_par_mape
            model_pars_df$det_par_mape <- det_par_mape
            model_pars_df$neg_log_like <- test.formula@negLogLike

  
            write.csv(model_pars_df, paste0("results/simulation_experiments/raw/model_parameters/", spec_name,"_", method_name, ".csv"), row.names=FALSE)


        }


            
        occ_par_list <- model_list[[method_name]][["occu_parameters"]]
        det_par_list <- model_list[[method_name]][["det_parameters"]]


        # check occupancy probability at each occupied site
        test.df$occupied_prob_est <- calculate_weighted_sum(occ_par_list, test.df)
        test.df$occupied_prob_est <- expit(test.df$occupied_prob_est)
        
        test.df$det_prob_est <- calculate_weighted_sum(det_par_list, test.df)
        test.df$det_prob_est <- expit(test.df$det_prob_est)
        
        pred_observ <- unlist(test.df$occupied_prob_est*test.df$det_prob_est)

        dets <- pred_observ[test.df$species_observed == TRUE]
        nondets <- pred_observ[test.df$species_observed == FALSE]

        # ROC Curve    
        roc <- roc.curve(scores.class0 = dets, scores.class1 = nondets, curve = F)

        # PR Curve
        pr <- pr.curve(scores.class0 = dets, scores.class1 = nondets, curve = F)
        
        occ.mape.i <- mean(abs((test.df$occupied_prob - test.df$occupied_prob_est)/test.df$occupied_prob)) * 100
        det.mape.i <- mean(abs((test.df$det_prob - test.df$det_prob_est)/test.df$det_prob)) * 100


        
        predictions_df <- data.frame(
            checklist_id = test.df$checklist_id,
            species_observed = test.df$species_observed,
            prediction = pred_observ
        )
        write.csv(predictions_df, paste0("results/simulation_experiments/raw/predictions/", spec_name,"_", method_name, "_run=", exp_run,".csv"), row.names=FALSE)

        

        metrics_df <- data.frame(
            run = exp_run,
            species = spec_name,
            method = method_name,
            occ_prob_mape = occ.mape.i,
            det_prob_mape = det.mape.i,
            auc = roc$auc,
            auprc = pr$auc.integral
        )
        
        write.csv(metrics_df, paste0("results/simulation_experiments/raw/metrics/", spec_name,"_", method_name, "_run=", exp_run,".csv"), row.names=FALSE)

    

        cat(paste0(method_name, " has AUC: ", round(roc$auc, 4), " and AUPRC: ", round(pr$auc.integral, 4), "\n"))


    }

}
