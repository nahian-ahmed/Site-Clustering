#############################################
# Summarize semi-simulated species results

# December 4, 2024
#############################################

library(dplyr)

# Initialize the summarized results dataframes


best_pars_summarized <- data.frame(matrix(ncol = 0, nrow = 0))

results_clust <- data.frame(matrix(ncol = 0, nrow = 0))

results_sum <- data.frame(matrix(ncol = 0, nrow = 0))

data_loss <- data.frame(matrix(ncol = 0, nrow = 0))

placeholder_spec_name <- "Sim-1"


for (spec_name in species_names) {

    for (m_name in method_names){
        
        method_name <- parseMethodName(m_name)
        


        if (startsWith(method_name, "BayesOpt")){
            
            bo_pars <- read.delim(paste0("results/simulation_experiments/raw/clustering_parameters/", spec_name,"_", method_name,".csv"), sep = ",", header = T)
            n_pars <- nrow(bo_pars)

            t_df <- data.frame(
                species = rep(spec_name, n_pars),
                method = rep(method_name, n_pars),
                par_name = bo_pars$parameter,
                best_par = bo_pars$value
            )
            best_pars_summarized <- rbind(best_pars_summarized, t_df)

        }

        clust_metrics <- read.delim(paste0("results/simulation_experiments/raw/clusterings/", spec_name,"_", method_name,"_clust_stats.csv"), sep = ",", header = T)
        

        t_df <- data.frame(
                species = spec_name,
                ref_method = species_df[species_df$species==spec_name,]$ref_method,
                method = method_name,
                ARI = clust_metrics$ari,
                AMI = clust_metrics$ami,
                NID = clust_metrics$nid
            )
        results_clust <- rbind(results_clust, t_df)


        
        for (exp_run in seq(1:runs)){

            perf_metrics <- read.delim(paste0("results/simulation_experiments/raw/metrics/", spec_name,"_", method_name, "_run=", exp_run,".csv"), sep = ",", header = T)
            
            model_pars <- read.delim(paste0("results/simulation_experiments/raw/model_parameters/", spec_name,"_", method_name, ".csv"), sep = ",", header = T)

            t_df <- data.frame(
                species = spec_name,
                ref_method = species_df[species_df$species==spec_name,]$ref_method,
                occ_rate = species_df_formatted[species_df_formatted$species==spec_name,]$occ_rate,
                det_rate =  species_df_formatted[species_df_formatted$species==spec_name,]$det_rate,
                prevalence = species_df_formatted[species_df_formatted$species==spec_name,]$prevalence,
                method = method_name,
                run = exp_run,
                ARI = clust_metrics$ari,
                AMI = clust_metrics$ami,
                NID = clust_metrics$nid
             )

            t_df$occ_intercept <- model_pars$occ_intercept
            for (occ_cov in occ_covs){
                t_df[, occ_cov] <- model_pars[, occ_cov] 
            }
            t_df$det_intercept <- model_pars$det_intercept
            for (det_cov in det_covs){
                t_df[, det_cov] <- model_pars[, det_cov] 
            }
                
            t_df$neg_log_like = model_pars$neg_log_like
            t_df$occ_par_mape = model_pars$occ_par_mape
            t_df$det_par_mape = model_pars$det_par_mape
            t_df$occ_prob_mape = perf_metrics$occ_prob_mape
            t_df$det_prob_mape = perf_metrics$det_prob_mape
            t_df$auc = perf_metrics$auc
            t_df$auprc = perf_metrics$auprc
            
            results_sum <- rbind(results_sum, t_df)

        

        }

        # Summarize data loss stats

        train.df.og <- read.delim(paste0("results/simulation_experiments/raw/data/", spec_name, "_train.csv"), sep = ",", header = T)
        og_train_size <- nrow(train.df.og)

        clustering <- read.delim(paste0("results/simulation_experiments/raw/clusterings/", spec_name,"_", method_name,".csv"), sep = ",", header = T)
        clust_on <- nrow(clustering)

        if ((clust_on < og_train_size) && (spec_name == placeholder_spec_name)){
            t_df <- data.frame(
                method = method_name,
                original_train_size = og_train_size,
                clustered_on = clust_on,
                thrown = og_train_size - clust_on,
                percentage_thrown = round((((og_train_size-clust_on)/og_train_size) * 100), 4)
                
            )
            data_loss <- rbind(data_loss, t_df)
        }

    }

    
}

# results_sum <- results_sum[order(results_sum$species, results_sum$method),]

# Write the summarized results to CSV
write.csv(best_pars_summarized, "results/simulation_experiments/summarized/best_pars_summarized.csv", row.names=FALSE)
write.csv(results_clust, "results/simulation_experiments/summarized/results_clust_summarized.csv", row.names=FALSE)
write.csv(results_sum, "results/simulation_experiments/summarized/results_summarized.csv", row.names=FALSE)
write.csv(data_loss, "results/simulation_experiments/summarized/data_loss_summarized.csv", row.names=FALSE)

