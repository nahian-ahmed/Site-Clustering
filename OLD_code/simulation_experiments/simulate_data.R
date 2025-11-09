########################################################
# Simulate species given reference clustering method

# December 4, 2024
########################################################




library(stringr)
library(ggplot2)
library(viridis)
library(ggpubr) # ggarrange function
library(ggnewscale) # for multiple scales
library(ggforce) # drawing ellipses over clusters
library(patchwork) # arranging plots of multiple sizes on same page
library(terra) # raster operations

dirs_to_create = c(
    "results",
    "results/simulation_experiments",
    "results/simulation_experiments/raw",
    "results/simulation_experiments/raw/data",
    "results/simulation_experiments/raw/test_data_splits",
    "results/simulation_experiments/raw/clusterings",
    "results/simulation_experiments/raw/clustering_parameters",
    "results/simulation_experiments/raw/model_parameters",
    "results/simulation_experiments/raw/predictions",
    "results/simulation_experiments/raw/metrics",
    "results/simulation_experiments/summarized",
    "results/simulation_experiments/plotted",
    "results/simulation_experiments/plotted/plots",
    "results/simulation_experiments/plotted/maps"
)

for (f.path in dirs_to_create){
    if (!dir.exists(f.path)) {
        dir.create(f.path)
    }
}






OR.init <- terra::rast("occupancy_feature_raster/occupancy_features.tif")
names(OR.init) <- occ_covs_all
OR.init <- subset(OR.init, occ_covs)


OR.normalized <- terra::rast("occupancy_feature_raster/occupancy_features_normalized.tif")
names(OR.normalized) <- occ_covs_all
OR.normalized <- subset(OR.normalized, occ_covs)



parse_method_names <- function(df){

    
    df$ref_method_parsed <- df$ref_method

    # Parse kmSq method names
    
    kmSq_idx <- str_detect(df$ref_method, 'kmSq')

    
    df[kmSq_idx, c('clustering_parameter', 'ref_method_parsed')] <- str_split_fixed(df[kmSq_idx,]$ref_method, '-', 2)

    # Calculate side length of each cell in km from area of cell in sqKm, and then convert from km to m    
    df$kmsq_par <- as.integer(sqrt(as.numeric(as.character(df$clustering_parameter))) * 1000)

    df$side_length_km <- sqrt(as.numeric(as.character(df$clustering_parameter)))

    df[kmSq_idx, ]$ref_method_parsed <- paste(df[kmSq_idx, ]$ref_method_parsed, "-", df[kmSq_idx, ]$kmsq_par, sep="")

    df$side_length_km <- round(df$side_length_km, 4)
    df[kmSq_idx, ]$clustering_parameter <- paste("kmSq=",df[kmSq_idx, ]$clustering_parameter, sep="")


    # Parse clustGeo method names
    clustGeo_idx <- str_detect(df$ref_method, 'clustGeo')


    df[clustGeo_idx,]$ref_method_parsed <- str_replace_all(df[clustGeo_idx, ]$ref_method_parsed, '-', '_')

    df[clustGeo_idx, c('ref_method_parsed', 'alpha', 'lambda')] <- str_split_fixed(df[clustGeo_idx,]$ref_method, '-', 3)


    df[clustGeo_idx, ]$ref_method_parsed <- paste(df[clustGeo_idx, ]$ref_method_parsed, "_", df[clustGeo_idx, ]$alpha, "_", df[clustGeo_idx, ]$lambda, sep="")

    df$alpha <- as.numeric(as.character(df$alpha)) / 100
    df$lambda <- as.numeric(as.character(df$lambda)) / 100 

    df[clustGeo_idx, ]$clustering_parameter <- paste("α=", df[clustGeo_idx, ]$alpha, ", λ=", df[clustGeo_idx, ]$lambda,sep="")

    df$clustering_parameter[is.na(df$clustering_parameter)] <- "-"
    df$side_length_km[is.na(df$side_length_km)] <- "-"

    df <- subset( df, select = -c( kmsq_par, alpha, lambda) )

    return (df)

}

prep_train_data <- function (placeholder_spec_name = "AMCR"){

    f.name_train <- paste0(placeholder_spec_name, "/", placeholder_spec_name, "_zf_filtered_region_2017.csv")

    train.df.og <- read.delim(paste0("checklist_data/species/", f.name_train), sep = ",", header = T)
  
    train.df.og <- train.df.og[!is.na(train.df.og$duration_minutes),]

    train.df.og <- train.df.og[train.df.og$observation_date >= "2017-05-15" & train.df.og$observation_date <= "2017-07-09",]
   
    train.df <- train.df.og


    train_env.df <- extractEnvFeat(train.df, OR.init, obs_covs = det_covs)
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



gen_train_data <-  function (sites_df, occ_par_list, det_par_list){
    
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

gen_test_data <- function (norm.list, occ_par_list, det_par_list, placeholder_spec_name = "AMCR"){
    
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



plot_gt_occu_maps <- function (species_df){



  
    plots <- list()
    i <- 1

    species_df$occ_prob_mean <- -1
    species_df$occ_prob_sd <- -1

    for (i in seq_len(nrow(species_df))){
    
        
        s_name = species_df[i,]$species
        m_name = species_df[i,]$ref_method


        par_lists <- get_parameters(species_df, i, occ_covs, det_covs)


        occ_par_list <- par_lists$occ_par_list


        predicted_map <- predict_sdm_map(occ_par_list, OR.normalized)

        predicted_occu_prob <- predicted_map$occ_prob


        names(predicted_occu_prob) <- c("occu")

        predicted_occu_prob <- aggregate(predicted_occu_prob, fact=6)

        predicted_occu_prob <- as.data.frame(predicted_occu_prob, xy=TRUE)
        

        species_df[i,]$occ_prob_mean <- round(mean(predicted_occu_prob$occu, na.rm = TRUE), 4)
        species_df[i,]$occ_prob_sd <- round(sd(predicted_occu_prob$occu, na.rm = TRUE), 4)

        x_min = min(predicted_occu_prob$x)
        x_max = max(predicted_occu_prob$x)
        y_min = min(predicted_occu_prob$y)
        y_max = max(predicted_occu_prob$y)

    
        pl <- ggplot() +
            geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgray", show.legend = FALSE)  +
            geom_raster(data = predicted_occu_prob, aes(x = x, y = y, fill = `occu`)) +
            scale_fill_viridis_c(option = "B", limits = c(0.0, 1.0)) +
            theme_void() +
            coord_fixed() +
            theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
            theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt")) +
            theme(
                legend.text = element_text(size = 14), 
                legend.title = element_text(size = 15, margin = margin(r=10, b=20)), 
                legend.key.size = unit(1, 'cm')
                ) +
            labs(fill = "Occupancy Probability")


        plots[[i]] <- pl

        i <- i + 1

        rm(pl)
        rm(predicted_occu_prob)
        
    }
    subplot_labels <- paste0(species_df$species,"; ", species_df$ref_method, "; Mean:",round(species_df$occ_prob_mean, 2),"; Prev.:",species_df$prevalence*100,"%")
    label_x_pos <- c(-0.3,-0.32,-0.345,-0.35,-0.425,-0.425,-0.45,-0.375,-0.275)
    
    nrow_c <- 3
    ncol_c <- 3
    plot_main <- ggarrange(
        plotlist = plots, 
        nrow = nrow_c, 
        ncol = ncol_c, 
        common.legend = TRUE, 
        legend = "bottom", 
        labels = subplot_labels, 
        font.label = list(size=12), 
        label.x = label_x_pos
    )
    
    
    height_p <- 18
    width_p <- 12
    
    png("results/simulation_experiments/plotted/maps/reference_occu_maps.png", units="in", width=width_p, height=height_p, res=300)

    
    show(plot_main)
    dev.off()

    return (species_df)
        
}


plotClustGroup <- function (species_df){

    zoomed_A = list(
        longitude = c(-123.025, -122.992),
        latitude = c(44.085, 44.118)
    )

    obs_plot <- NA
    plots <- list()

    valid_boundary <- terra::vect("occupancy_feature_raster/boundary/boundary.shp")
    crs(valid_boundary) <- crs(OR.init)
    
    region <- terra::crop(OR.init$elevation, valid_boundary, mask = TRUE)

    # region <- aggregate(region, fact=6)

    base_rast <- as.data.frame(region, xy=TRUE)
    
    base_rast_min <- min(base_rast$elevation)
    base_rast_max <- max(base_rast$elevation)
    
    x_min = min(base_rast$x)
    x_max = max(base_rast$x)
    y_min = min(base_rast$y)
    y_max = max(base_rast$y)

    i <- 1
    for (i in seq_len(nrow(species_df))){

        s_name = species_df[i,]$species

        pts_df <- read.delim(paste0("results/simulation_experiments/raw/data/", s_name,"_train.csv"), sep = ",", header = T)
        pts_df$site <- as.factor(pts_df$site)

        if (i==1){

            obs_plot <- ggplot()+ 
            
                geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgray", show.legend = FALSE)  +
                geom_raster(data = base_rast, aes(x = x, y = y, fill = `elevation`), show.legend = TRUE) +
                scale_fill_viridis_c(option = "H", limits = c(base_rast_min, base_rast_max)) +
                geom_point(data = pts_df, aes(x=longitude, y=latitude), color="black", size= 1, shape= 21, fill = "#FBFAF5", show.legend = FALSE) +
                annotate(
                    "segment", 
                    x = -123.75, 
                    xend = zoomed_A$longitude[1], 
                    y = 44.425, 
                    yend = (zoomed_A$latitude[1] + zoomed_A$latitude[2])/2, 
                    colour = "yellow"
                    ) +
                annotate(
                    "label", 
                    x = -124.05, 
                    y = 44.425, 
                    label = "Zoomed-in region",
                    size = 3.25,
                    ) +
                geom_rect(aes(xmin = zoomed_A$longitude[1], xmax = zoomed_A$longitude[2], ymin = zoomed_A$latitude[1], ymax = zoomed_A$latitude[2]), fill = "red", color="red", alpha=0.1, show.legend = FALSE)  +
                labs(
                    title = "Species\nObservations",
                    x = "Longitude",
                    y = "Latitude",
                    fill = "Elevation (m)"
                ) +

                theme_bw() +
                coord_fixed() +
                theme(
                    plot.title = element_text(hjust = 0.5, vjust = -52.5, face="bold", size = 12)
                ) +
                theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) +
                theme(
                    legend.position = "inside",
                    legend.position.inside = c(0.5, -0.22),
                    legend.direction="horizontal",
                    legend.title = element_text(size = 12, margin = margin(r=5, b=20)), 
                    legend.key.size = unit(0.5, 'cm'),
                    legend.text = element_text(size = 8, angle = 90)
                    ) 
                


        }

        pts_df <- pts_df[
            (pts_df$longitude > zoomed_A$longitude[1]) &
            (pts_df$longitude < zoomed_A$longitude[2]) &
            (pts_df$latitude > zoomed_A$latitude[1]) &
            (pts_df$latitude < zoomed_A$latitude[2]) ,]

        pts_df$lat_long <- paste0(pts_df$latitude, "_", pts_df$longitude)
        pts_df <- dplyr::distinct(pts_df, lat_long, .keep_all = T)
        
        plots[[i]] <- ggplot() +
            geom_rect(
                aes(
                    xmin = zoomed_A$longitude[1], 
                    xmax = zoomed_A$longitude[2], 
                    ymin = zoomed_A$latitude[1], 
                    ymax = zoomed_A$latitude[2]
                ), 
                fill = "#E6E6E6", 
                color = "#E6E6E6", 
                show.legend = FALSE
            ) +
            geom_raster(
                data = base_rast, 
                aes(x = x, y = y, fill = elevation),
                show.legend = FALSE
            ) +
            scale_fill_viridis_c(
                option = "H",
                aesthetics = "fill"
            ) +
            new_scale_fill() + 
            geom_point(
                data = pts_df, 
                aes(
                    x = longitude, 
                    y = latitude,
                    fill = site 
                ), 
                shape = 21,
                size = 3.5, 
                color = "black",
                show.legend = FALSE
            ) +
            geom_mark_ellipse(
                data = pts_df,
                aes(
                    x = longitude,
                    y = latitude, 
                    fill = site, 
                    color = site 
                ),
                alpha = 0.1,
                show.legend = FALSE
            ) +
            scale_fill_discrete() +
            scale_color_discrete() +
            coord_fixed(
                xlim = c(zoomed_A$longitude[1], zoomed_A$longitude[2]),
                ylim = c(zoomed_A$latitude[1], zoomed_A$latitude[2])
            ) +
            theme_bw() +
            labs(
                x = "Longitude",
                y = "Latitude"
            ) +
            theme(plot.margin = margin(t = 0.25, r = 0.0, b = 0.1, l = 0.0, "in"))

        


        i <- i + 1

        
    }
    
 
    label_x_pos <- c( 0.2, 0.17, 0.125, 0.1, 0.01, 0.01, 0.01, 0.01, 0.25)
    species_df$labels <- paste0(species_df$species, "; ", species_df$ref_method)    

    plot_clust <- ggarrange(
        plotlist=plots, 
        nrow=3, 
        ncol=3, 
        labels=species_df$labels,
        label.x=label_x_pos,
        font.label = list(size = 12, face = "bold", vjust = -1),
        common.legend = TRUE
         
    )

    plot_main <- 
        obs_plot + 
        plot_clust  + 
        plot_layout(nrow = 1, byrow = TRUE, widths = c(1, 4.25))

 
    height_p <- 8
    width_p <- 12

    png(paste0("results/simulation_experiments/plotted/maps/reference_clust_maps.png"), units="in", width=width_p, height=height_p, res=300)
    show(plot_main)
    dev.off()


        
}





simulate_species <- function(df){


    df <- parse_method_names(df)
    stat_colnames <- c(
        "occ_rate", 
        "det_rate", 
        "prevalence", 
        "num_clusters", 
        "min_cluster_size", 
        "max_cluster_size", 
        "mean_cluster_size", 
        "sd_cluster_size", 
        "percentage_svs", 
        "avg_intra_site_var"
    )
    df[stat_colnames] <- NA

    train.df.sup <- prep_train_data()
    train.df <- train.df.sup$train.df
    norm.list <- train.df.sup$norm.list

    n_species <- nrow(df)

    for (i in seq_len(n_species)){
        

        set.seed(1)

        method_names <- list(df[i,]$ref_method_parsed)


        species_name <- df[i,]$species


        methods <- genExp(method_names)
        groupedSite <- getClusterings(methods, train.df, occ_covs, det_covs)

        
        for(method_name in names(groupedSite)){


            groupedSite_df <- groupedSite[[method_name]]
            sites_occ <- subset(groupedSite_df, !duplicated(site))$site
            closed_df <- enforceClosure(groupedSite_df, occ_covs, sites_occ)

            par_lists <- get_parameters(df, i, occ_covs, det_covs)
    
            occ_par_list <- par_lists$occ_par_list
            det_par_list <- par_lists$det_par_list

            species_train_df <- gen_train_data(closed_df, occ_par_list, det_par_list)
            species_test_df <- gen_test_data(norm.list, occ_par_list, det_par_list)

            # Write normalized train and test data to csv
            write.csv(species_train_df, paste0("results/simulation_experiments/raw/data/", species_name,"_train.csv"))
            write.csv(species_test_df, paste0("results/simulation_experiments/raw/data/", species_name,"_test.csv"))

            
            # Calculate clustering and species specific stats
            clust_descr_stats <- calcDescriptiveClusteringStatsWithReference(species_train_df, "site", occ_covs, normalize = TRUE)

            df[i,]$num_points <- clust_descr_stats$n_points
            df[i,]$num_clusters <- clust_descr_stats$n_clusters 
            df[i,]$min_cluster_size <- clust_descr_stats$min_size
            df[i,]$max_cluster_size <- clust_descr_stats$max_size
            df[i,]$mean_cluster_size <- clust_descr_stats$mean_size
            df[i,]$sd_cluster_size <- clust_descr_stats$sd_size
            df[i,]$percentage_svs <- clust_descr_stats$perc_svs
            df[i,]$avg_intra_site_var <- clust_descr_stats$average_intra_cluster_variance
            
            df[i,]$occ_rate <- round(sum(species_train_df$occupied)/nrow(species_train_df), 4)
            df[i,]$det_rate <- round(sum(species_train_df$detection)/nrow(species_train_df), 4)
            df[i,]$prevalence <- round(sum(species_train_df$species_observed)/nrow(species_train_df), 4)
            
            
            
        }


    }


    # plot maps of all simulated species
    df <- plot_gt_occu_maps(df)
    
    plotClustGroup(df)

    # remove unnecessary columns
    df <- subset(df, select = -c(ref_method_parsed))
    
    # write summaries of species
    write.csv(df, "results/simulation_experiments/summarized/simulation_species_formatted.csv")


    
}

simulate_species(species_df)