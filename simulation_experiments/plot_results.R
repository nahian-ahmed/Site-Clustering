#########################
# Plot results and maps

# Octcober 8, 2024
#########################

library(reshape2)
library(viridis)
library(ggplot2)
library(ggrepel)
library(ggpubr) # ggarrange function
library(ggh4x) # for panel size configuration
library(ggnewscale) # for multiple scales
library(ggforce) # drawing ellipses over clusters
library(cowplot) # for placing grobs next to each other
library(patchwork) # arranging plots of multiple sizes on same page
library(lmerTest) # linear mixed effect models
library(terra) # raster operations
library(dplyr) # For data manipulation
library(geosphere) # For geographic calculations (distance and area)

plot_maps <- FALSE


OR.init <- terra::rast("occupancy_feature_raster/occupancy_features.tif")
names(OR.init) <- occ_covs_all
OR.init <- subset(OR.init, occ_covs)


OR.normalized <- terra::rast("occupancy_feature_raster/occupancy_features_normalized.tif")
names(OR.normalized) <- occ_covs_all
OR.normalized <- subset(OR.normalized, occ_covs)


species_descr_ext <- read.delim("results/simulation_experiments/summarized/simulation_species_formatted.csv", sep = ",", header = T)
species_descr_ext$species_ext <- paste0(species_descr_ext$ref_method, "\nOcc. Rate: ", species_descr_ext$occ_rate, "; Prev.: ", species_descr_ext$prevalence)


results_sum_a <- read.delim("results/simulation_experiments/summarized/results_summarized.csv", sep = ",", header = T)

dropped_colnames <- c("occ_intercept", occ_covs, "det_intercept", det_covs, "neg_log_like")

results_sum_a <- results_sum_a[,!(colnames(results_sum_a) %in% dropped_colnames)]

colors <- c("red", "navy", "cyan", "pink", "green", "brown", "purple", "yellow", "blue", "darkgrey", "forestgreen")

shapes <- c(0, 1, 2, 12, 10, 5, 6, 7, 9)

species_all <- unique(results_sum_a$species)

alg_all <- unique(results_sum_a$method)



alg_all_o <- c(
    "reference-clustering",
    "1to10",
    "2to10",
    "2to10-sameObs",
    "1-kmSq",
    "lat-long",
    "rounded-4",
    "SVS",
    "1-UL",
    "DBSC",
    "BayesOptClustGeo"
)



alg_all_o_complete <- c(
    "reference-clustering",
    "1-kmSq",
    "lat-long",
    "rounded-4",
    "DBSC",
    "BayesOptClustGeo"
)

alg_all_o_mv <- c(
    "reference-clustering",
    "1to10",
    "2to10",
    "2to10-sameObs",
    "1-kmSq",
    "lat-long",
    "rounded-4",
    "DBSC",
    "BayesOptClustGeo"
)

alg_all_o_cl_hist <- c(
    "1-kmSq",
    "rounded-4",
    "DBSC",
    "BayesOptClustGeo"
)

metrics_list_formatted <- list(
    ARI = "Adjusted Rand Index",
    AMI = "Adjusted Mutual Information",
    NID = "Normalized Information Difference",
    occ_rate = "Occupancy Rate",
    det_rate = "Detection Rate",
    prevalence = "Prevalence",
    occ_par_mape = "Occupancy Parameter MAPE",
    det_par_mape = "Detection Parameter MAPE",
    occ_prob_mape = "Occupancy Probability MAPE",
    det_prob_mape = "Detection Probability MAPE",
    auc = "AUC",
    auprc = "AUPRC"
)

repeat_nonvarying_metrics = c("ARI", "AMI", "NID", "occ_par_mape",  "det_par_mape")
repeat_varying_metrics = c("occ_prob_mape", "det_prob_mape", "auc", "auprc")
clustering_metrics = c("ARI", "AMI", "NID")
species_characteristics_metrics = c("occ_rate", "det_rate", "prevalence")
primary_performance_metrics = c("occ_par_mape", "det_par_mape")
metric_names <- names(metrics_list_formatted)

# 11:12 are the prediction metric columns
metric_names_downstream <- metric_names[9:12] 

for (metric_name in metric_names){
    
    names(results_sum_a)[names(results_sum_a) == metric_name] <- metrics_list_formatted[[metric_name]]
    
}



results_sum_a$species_ref_method <- paste0(results_sum_a$species,"\n(",results_sum_a$ref_method,")")



##############################################################
# 1. SDM maps
##############################################################


########################################
# Function to plot occupancy maps
########################################
plotMapGroup <- function (species, sp_i, alg_all){

    nrow_c <- 3
    height_p <- 24
    width_p <- 20
    
    png(paste0("results/simulation_experiments/plotted/maps/", species, "_occu_maps.png"), units="in", width=width_p, height=height_p, res=300)

    plots <- list()
    means <- list()
    # morans <- list()
    
    i <- 1
            
    par_lists <- get_parameters(species_df, sp_i, occ_covs, det_covs)

    occ_par_list <- par_lists$occ_par_list


    predicted_map <- predict_sdm_map(occ_par_list, OR.normalized)
    
    predicted_occu_prob <- predicted_map$occ_prob

    # morans[[i]] <- round(Moran(predicted_occu_prob), 4) 
    
    names(predicted_occu_prob) <- c("occu")

    means[[i]] <- round(mean(as.matrix(predicted_occu_prob$occu), na.rm = TRUE), 4)

    predicted_occu_prob <- aggregate(predicted_occu_prob, fact=6)

    predicted_occu_prob <- as.data.frame(predicted_occu_prob, xy=TRUE)
    
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
    
    
    for (m in alg_all){
        
        model_pars <- read.delim(paste0("results/simulation_experiments/raw/model_parameters/", species,"_", m, ".csv"), sep = ",", header = T)
            
        par_lists <- get_parameters(model_pars, 1, occ_covs, det_covs)

        occ_par_list <- par_lists$occ_par_list


        predicted_map <- predict_sdm_map(occ_par_list, OR.normalized)
        
        predicted_occu_prob <- predicted_map$occ_prob

        # morans[[i]] <- round(Moran(predicted_occu_prob), 4) 
        
        names(predicted_occu_prob) <- c("occu")

        means[[i]] <- round(mean(as.matrix(predicted_occu_prob$occu), na.rm = TRUE), 4)

        predicted_occu_prob <- aggregate(predicted_occu_prob, fact=6)

        predicted_occu_prob <- as.data.frame(predicted_occu_prob, xy=TRUE)




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
    
    
    # alg_all[alg_all=="ground-truth"] <- "ground-truth-clustering"

    subplot_labels <- paste0(c("reference-occupancy",alg_all),"; Mean: ",unlist(means))
    label_x_pos <- c(-0.225, -0.2475, 0.1, 0.1, -0.125, 0.075, 0.05, 0.0, 0.125, 0.125, 0.1, -0.2)
    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=4, common.legend=TRUE, legend="bottom", labels=subplot_labels, font.label=list(size=16), label.x=label_x_pos)
        
    
    show(plot_main)
    dev.off()

        
}



########################################
# Function to plot occupancy difference maps
########################################
plotDiffMapGroup <- function (species, sp_i, alg_all){

    nrow_c <- 3
    ncol_c <- 4
    height_p <- 24
    width_p <- 20
    
    png(paste0("results/simulation_experiments/plotted/maps/", species, "_occu_diff_maps.png"), units="in", width=width_p, height=height_p, res=300)

    plots <- list()
    
    i <- 1
            
    par_lists <- get_parameters(species_df, sp_i, occ_covs, det_covs)
    occ_par_list <- par_lists$occ_par_list

    predicted_map <- predict_sdm_map(occ_par_list, OR.normalized)
    predicted_occu_prob <- predicted_map$occ_prob
    names(predicted_occu_prob) <- c("occu")
    
    predicted_occu_prob_gt <- predicted_occu_prob
    predicted_occu_prob <- aggregate(predicted_occu_prob, fact=6)

    predicted_occu_prob <- as.data.frame(predicted_occu_prob, xy=TRUE)
    
    x_min <- min(predicted_occu_prob$x)
    x_max <- max(predicted_occu_prob$x)
    y_min <- min(predicted_occu_prob$y)
    y_max <- max(predicted_occu_prob$y)
    
    # First plot (True Occupancy)
    pl <- ggplot() +
            geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgray", show.legend = FALSE)  +
            geom_raster(data = predicted_occu_prob, aes(x = x, y = y, fill = `occu`)) +
            scale_fill_viridis_c(option = "B", limits = c(0.0, 1.0), name = "Occupancy Probability") +
            theme_void() +
            coord_fixed() +
            theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
            theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
            theme(legend.text = element_text(size = 14), legend.title = element_text(size = 15, margin = margin(r=10, b=20)), legend.key.size = unit(1, 'cm'), legend.direction = "horizontal") +
            labs(fill = "Occupancy Probability")

    plots[[i]] <- pl
    i <- i + 1

    legend1 <- get_legend(pl)  # Extract legend from first plot

    for (m in alg_all){
        
        model_pars <- read.delim(paste0("results/simulation_experiments/raw/model_parameters/", species,"_", m, ".csv"), sep = ",", header = T)
            
        par_lists <- get_parameters(model_pars, 1, occ_covs, det_covs)
        occ_par_list <- par_lists$occ_par_list

        predicted_map <- predict_sdm_map(occ_par_list, OR.normalized)
        predicted_occu_prob <- predicted_map$occ_prob
        names(predicted_occu_prob) <- c("occu")

        predicted_occu_prob <- get_occu_map_diff(predicted_occu_prob_gt, predicted_occu_prob)
        
        diff_df <- data.frame(
            species = species,
            method = m,
            perc_diff_mean = round(mean(as.matrix(predicted_occu_prob$percentage_difference), na.rm = TRUE), 6),
            perc_diff_sd = round(sd(as.matrix(predicted_occu_prob$percentage_difference), na.rm = TRUE), 6)
        )
        
        write.csv(diff_df, paste0("results/simulation_experiments/raw/metrics/occu_map_perc_diff_", species,"_", m,".csv"), row.names=FALSE)

        predicted_occu_prob <- aggregate(predicted_occu_prob, fact=6) 

        predicted_occu_prob <- as.data.frame(predicted_occu_prob, xy=TRUE)

        pl <- ggplot() +
            geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgray", show.legend = FALSE)  +
            geom_raster(data = predicted_occu_prob, aes(x = x, y = y, fill = `difference`)) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1.0, 1.0), name = "Occupancy Probability Difference") +
            theme_void() +
            coord_fixed() +
            theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
            theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
            theme(legend.text = element_text(size = 14), legend.title = element_text(size = 15, margin = margin(r=10, b=20)), legend.key.size = unit(1, 'cm'), legend.direction = "horizontal") +
            labs(fill = "Occupancy Probability Difference")

        plots[[i]] <- pl


        if (i == 2) {
            legend2 <- get_legend(pl)  # Extract legend from second plot
        }
        i <- i + 1
    }
    
    # alg_all[alg_all=="ground-truth"] <- "ground-truth-clustering"
    subplot_labels <- c("reference-occupancy", alg_all)

    label_x_pos <- c(0.1, 0.135, 0.4, 0.4, 0.2, 0.35, 0.375, 0.325, 0.45, 0.435, 0.435, 0.125)

    # Remove legends from individual plots
    plots <- lapply(plots, function(p) p + theme(legend.position="none"))

    # Arrange plots without legends
    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=ncol_c, labels=subplot_labels, font.label=list(size=16), label.x=label_x_pos)

    # Arrange legends horizontally and place them side by side
    legend_combined <- plot_grid(legend1, legend2, ncol=2, rel_widths=c(1, 2), align="h")

    # Combine everything
    final_plot <- plot_grid(plot_main, legend_combined, ncol=1, rel_heights=c(5, 0.2))

    show(final_plot)
    dev.off()
}


########################################
# Function to plot clustered sites
########################################

plotClustGroup <- function (species, alg_all){
    

    zoomed_A = list(
        longitude = c(-123.025, -122.992),
        latitude = c(44.085, 44.118)
    )
    obs_plot <- NA
    plots <- list()

    valid_boundary <- terra::vect("occupancy_feature_raster/boundary/boundary.shp")
    crs(valid_boundary) <- crs(OR.normalized)
    
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
    for (m in alg_all){

        pts_df <- read.delim(paste0("results/simulation_experiments/raw/clusterings/", species,"_", m, ".csv"), sep = ",", header = T)
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
                    plot.title = element_text(hjust = 0.5, vjust = -25, face="bold", size = 12)
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
 
    label_x_pos <- c(0.075, 0.4, 0.34, 0.45, 0.1)
    

    plot_clust <- ggarrange(
        plotlist=plots, 
        nrow=2, 
        ncol=3, 
        labels=alg_all,
        label.x=label_x_pos,
        font.label = list(size = 12, face = "bold", vjust = -1)
         
    )

    plot_main <- 
        obs_plot + 
        plot_clust  + 
        plot_layout(nrow = 1, byrow = TRUE, widths = c(1, 4.25))

 

    height_p <- 6
    width_p <- 12

    png(paste0("results/simulation_experiments/plotted/maps/",species,"_clust_maps.png"), units="in", width=width_p, height=height_p, res=300)
    show(plot_main)
    dev.off()


        
}

#  # Plot all grouped maps

# species_all <- c("Sim-1")

alg_clust <- c("reference-clustering", "1-kmSq", "rounded-4", "DBSC", "BayesOptClustGeo")

if (plot_maps){
    perc_diff_df <- data.frame(matrix(nrow=0, ncol=0))

    sp_i <- 1
    for (s in species_all){

        #######
        # 1.1. Plot Occupancy Maps
        #######  
        plotMapGroup(s, sp_i, alg_all_o)
        
        #######
        # 1.2. Plot Occupancy Difference Maps
        #######
        plotDiffMapGroup(s, sp_i, alg_all_o)

        #######
        # 1.3. Plot Cluster Maps
        #######
        plotClustGroup(s, alg_clust)

        for (m in alg_all_o){
            t_df <- read.delim(paste0("results/simulation_experiments/raw/metrics/occu_map_perc_diff_", s,"_", m,".csv"), sep = ",", header = T)
            perc_diff_df <- rbind(perc_diff_df, t_df)
        }
        sp_i <- sp_i + 1
    }

    write.csv(perc_diff_df, paste0("results/simulation_experiments/summarized/occu_map_perc_diff_summarized.csv"), row.names=FALSE)

}



#######
# 1.4. Percentage Occu Map Diff Barplot
#######

height_p <- 9
width_p <- 9
nrow_c <- 3
ncol_c <- 3

perc_diff_df <- read.delim("results/simulation_experiments/summarized/occu_map_perc_diff_summarized.csv", sep = ",", header = T)

png(paste0("results/simulation_experiments/plotted/plots/metric=occ-map-mape.png"), units="in", width=width_p, height=height_p, res=300)


i <- 1
plots <- list()

# Iterate over species
for (spec in species_names){

    perc_diff_df_c <- perc_diff_df[perc_diff_df$species==spec,]
    perc_diff_df_c$method <- factor(perc_diff_df_c$method, levels = alg_all_o) 
        
            
    pl <- ggplot(perc_diff_df_c, aes(x=method, y=perc_diff_mean, fill = method)) +
        geom_bar(stat="identity",color = "black") +
        geom_errorbar(aes(ymin=(perc_diff_mean - perc_diff_sd), ymax=(perc_diff_mean + perc_diff_mean)), width=.2,
            position=position_dodge(.9)) +
        scale_fill_manual(name = "Algorithm", labels= alg_all_o, values=colors, drop=FALSE, na.translate=FALSE) +
        theme_linedraw() +
        labs(title = spec, subtitle = species_descr_ext[species_descr_ext$species==spec,]$species_ext) +
        theme(
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.title.y = element_blank(), 
            plot.title = element_text(size = 11, hjust=0.5, face = "bold"),
            plot.subtitle = element_text(size = 10, hjust = 0.5)
        )    
    
    plots[[i]] <- pl

    i <- i + 1

}

plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=ncol_c, common.legend=TRUE, legend="bottom")

plot_main <- annotate_figure(plot_main,
        left = text_grob("Occupancy Probability Map MAPE", size = 13, hjust = 0.3, color = "black", rot = 90),
)    


show(plot_main)
dev.off()




#######################################################################
# 2. Clustering Quality, Model Parameter, Prediction Error Plots
#######################################################################


width_p <- 9
height_p <- 9
nrow_c <- 3
ncol_c <- 3

# Iterate over performance metrics

for (metric_name in repeat_nonvarying_metrics){
    
    
    #############################
    # 2.1. All methods
    #############################

    png(paste0("results/simulation_experiments/plotted/plots/metric=", gsub("_", "-", metric_name),".png"), units="in", width=width_p, height=height_p, res=300)

    i <- 1
    plots <- list()

    # Iterate over species
    for (spec in species_names){
    
        results_sum <- results_sum_a[results_sum_a$species==spec & results_sum_a$run==1 ,]
        results_sum$method <- factor(results_sum$method, levels = alg_all_o) 
            
        pl <- ggplot(results_sum, aes(x=method, y=!!as.name(metrics_list_formatted[[metric_name]]), fill = method)) +
            geom_bar(stat="identity",color = "black") +
            scale_fill_manual(name = "Algorithm", labels= alg_all_o, values=colors, drop=FALSE, na.translate=FALSE) +
            theme_linedraw() +
            labs(title = spec, subtitle = species_descr_ext[species_descr_ext$species==spec,]$species_ext) +
            theme(
                axis.title.x = element_blank(), 
                axis.text.x = element_blank(), 
                axis.title.y = element_blank(), 
                plot.title = element_text(size = 11, hjust=0.5, face = "bold"),
                plot.subtitle = element_text(size = 10, hjust = 0.5)
            )
        
        
        plots[[i]] <- pl

        i <- i + 1

    }
    
    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=ncol_c, common.legend=TRUE, legend="bottom")

    plot_main <- annotate_figure(plot_main,
            left = text_grob(metrics_list_formatted[[metric_name]], size = 13, hjust = 0.3, color = "black", rot = 90),
    )    

    
    show(plot_main)
    dev.off()


  

}



results_sum_prediction <- results_sum_a

results_sum_prediction <- results_sum_prediction[, !(colnames(results_sum) %in% repeat_nonvarying_metrics)]



for (metric_name in metric_names_downstream){ 
    results_sum_prediction[,paste0(metrics_list_formatted[[metric_name]]," Mean")] <- NA

}

# Calculate mean and sd
for (s_name in species_names){
    for (alg in alg_all_o){
        
        c_df <- results_sum_prediction[((results_sum_prediction$species==s_name) & (results_sum_prediction$method==alg)),]

        
        for (metric_name in metric_names_downstream){
            
        
            results_sum_prediction[((results_sum_prediction$species==s_name) & (results_sum_prediction$method==alg)),paste0(metrics_list_formatted[[metric_name]]," Mean")] <- mean(c_df[,metrics_list_formatted[[metric_name]]], )
            results_sum_prediction[((results_sum_prediction$species==s_name) & (results_sum_prediction$method==alg)),paste0(metrics_list_formatted[[metric_name]]," SD")] <- sd(c_df[,metrics_list_formatted[[metric_name]]], )
        }

    }
}


for (metric_name in metric_names_downstream){
    results_sum_prediction <- results_sum_prediction[, !(colnames(results_sum_prediction) %in% c(metrics_list_formatted[[metric_name]]))]
    names(results_sum_prediction)[names(results_sum_prediction) == paste0(metrics_list_formatted[[metric_name]]," Mean")] <- metrics_list_formatted[[metric_name]]
    
}
for (metric_name in repeat_varying_metrics){
    
    #############################
    # 2.2. Downstream metrics
    #############################

    png(paste0("results/simulation_experiments/plotted/plots/metric=", gsub("_", "-", metric_name),".png"), units="in", width=width_p, height=height_p, res=300)

    i <- 1
    plots <- list()

    # Iterate over species
    for (spec in species_names){
        results_sum_prediction_c <- results_sum_prediction[results_sum_prediction$species==spec & results_sum_prediction$run==1 ,]
        results_sum_prediction_c$method <- factor(results_sum_prediction_c$method, levels = alg_all_o) 
            
        pl <- ggplot(results_sum_prediction_c, aes(x=method, y=!!as.name(metrics_list_formatted[[metric_name]]), fill = method)) +
            geom_bar(stat="identity",color = "black") +
            geom_errorbar(aes(ymin=!!as.name(metrics_list_formatted[[metric_name]]) - !!as.name(paste0(metrics_list_formatted[[metric_name]], " SD")), ymax=!!as.name(metrics_list_formatted[[metric_name]]) + !!as.name(paste0(metrics_list_formatted[[metric_name]], " SD"))), width=.2,
                position=position_dodge(.9)) +
            scale_fill_manual(name = "Algorithm", labels= alg_all_o, values=colors, drop=FALSE, na.translate=FALSE) +
            theme_linedraw() +
            labs(title = spec, subtitle = species_descr_ext[species_descr_ext$species==spec,]$species_ext) +
            theme(
                axis.title.x = element_blank(), 
                axis.text.x = element_blank(), 
                axis.title.y = element_blank(), 
                plot.title = element_text(size = 11, hjust=0.5, face = "bold"),
                plot.subtitle = element_text(size = 10, hjust = 0.5)
            )
            
        
        
        plots[[i]] <- pl

        i <- i + 1

    }
    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=ncol_c, common.legend=TRUE, legend="bottom")

    plot_main <- annotate_figure(plot_main,
            left = text_grob(metrics_list_formatted[[metric_name]], size = 13, hjust = 0.3, color = "black", rot = 90),
    )    


    show(plot_main)
    dev.off()

}





#########################################################################
# 3. Clustering Quality vs Model Parameter and Prediction Error Plots
#########################################################################

height_p <- 9
width_p <- 9.5
nrow_c <- 3
ncol_c <- 3

for (metric_name in clustering_metrics){



    alg_o <- alg_all_o_complete
    colors_c <- colors[!(colors %in% c("navy","blue","cyan","pink", "yellow"))]
    
    #############
    # 3.1. vs Occ Par MAPE
    #############

    png(paste0("results/simulation_experiments/plotted/plots/x=", gsub("_", "-", metric_name),"_y=occ-par-mape.png"), units="in", width=width_p, height=height_p, res=300)

    i <- 1
    plots <- list()

    results_sum_t <- results_sum_a[results_sum_a$method %in% alg_o,]

    sup_x_title <- metrics_list_formatted[[metric_name]]
    
    sup_x_title_hjust <- 0.75
    if(metric_name %in% c("NID")){
        results_sum_t[,metrics_list_formatted[[metric_name]]] <- -1 * results_sum_t[,metrics_list_formatted[[metric_name]]]
        sup_x_title <- paste0("Negative of ", sup_x_title, " (↑)")
        sup_x_title_hjust <- 0.6
    }
    else if (metric_name %in% c("ARI", "AMI")){
        sup_x_title <- paste0(sup_x_title, " (↑)")
    }

    results_sum_t[,"Occupancy Parameter MAPE"] <- -1 * results_sum_t[,"Occupancy Parameter MAPE"]

    # Iterate over species
    for (spec in species_names){
        results_sum <- results_sum_t[results_sum_t$species==spec & results_sum_t$run==1,]
        results_sum$method <- factor(results_sum$method, levels = alg_o) 

        x_min <- min(results_sum[,metrics_list_formatted[[metric_name]]])
        x_max <- max(results_sum[,metrics_list_formatted[[metric_name]]])
        y_min <- min(results_sum[,"Occupancy Parameter MAPE"])
        y_max <- max(results_sum[,"Occupancy Parameter MAPE"])
      
        slope_diag <- (y_max - y_min) / (x_max - x_min)
        int_diag <- y_max - (slope_diag * x_max)

      
        # Correlation test
        correlation_result <- cor.test(results_sum[,metrics_list_formatted[[metric_name]]], results_sum[,"Occupancy Parameter MAPE"], method = "pearson") 
        pearson_r <- format(round(correlation_result$estimate, 2), nsmall = 2)
        p_value <- format(round(correlation_result$p.value, 2), nsmall = 2)

      
        annotate_hjust <- -2.7
        if (pearson_r<0){
            annotate_hjust <- -2.45
        }

        pl <- ggplot(results_sum, aes(x=!!as.name(metrics_list_formatted[[metric_name]]), y=!!as.name("Occupancy Parameter MAPE"), color = method)) +
            geom_abline(
                slope = slope_diag,
                intercept = int_diag,
                color="#A9A9A9",
                linetype = "dashed"
            ) +
            geom_point(size=3, stroke=1, shape = 8) +
            theme_linedraw() +
            scale_color_manual(name = "Algorithm", labels = alg_o, values=colors_c, drop = FALSE, na.translate = FALSE) +
            labs(title = spec, subtitle = species_descr_ext[species_descr_ext$species==spec,]$species_ext) +
            theme(
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(), 
                plot.title = element_text(size = 11, hjust=0.5, face = "bold"),
                plot.subtitle = element_text(size = 10, hjust = 0.5)
            ) +
            annotate(geom="text", x=-Inf, y=-Inf, label=paste0("r = ", pearson_r),  hjust = annotate_hjust, vjust = -2, color="black") +
            annotate(geom="text", x=-Inf, y=-Inf, label=paste0("p = ", p_value),  hjust = -2.5, vjust = -0.45, color="black")
                
        
        
        plots[[i]] <- pl

        i <- i + 1

    }
    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=ncol_c, common.legend=TRUE, legend="right")

    plot_main <- annotate_figure(
        plot_main,
        left = text_grob("Negative of Occupancy Parameter MAPE (↑)", size = 13, hjust = 0.5, color = "black", rot = 90),
        bottom = text_grob(sup_x_title, size = 13,  hjust = sup_x_title_hjust, ,color = "black")
    )    


    show(plot_main)
    dev.off()

  

    #############
    # 3.2. vs Occ Prob MAPE
    #############

    png(paste0("results/simulation_experiments/plotted/plots/x=", gsub("_", "-", metric_name),"_y=occ-prob-mape.png"), units="in", width=width_p, height=height_p, res=300)

    i <- 1
    plots <- list()

    results_sum_prediction_t <- results_sum_prediction[results_sum_prediction$method %in% alg_o,]

    if(metric_name %in% c("NID")){
        results_sum_prediction_t[,metrics_list_formatted[[metric_name]]] <- -1 * results_sum_prediction_t[,metrics_list_formatted[[metric_name]]]

    }



    results_sum_prediction_t[,"Occupancy Probability MAPE"] <- -1 * results_sum_prediction_t[,"Occupancy Probability MAPE"]

    # Iterate over species
    for (spec in species_names){
        results_sum_prediction_c <- results_sum_prediction_t[results_sum_prediction_t$species==spec & results_sum_prediction_t$run==1,]
        results_sum_prediction_c$method <- factor(results_sum_prediction_c$method, levels = alg_o) 
        

        x_min <- min(results_sum_prediction_c[,metrics_list_formatted[[metric_name]]])
        x_max <- max(results_sum_prediction_c[,metrics_list_formatted[[metric_name]]])
        y_min <- min(results_sum_prediction_c[,"Occupancy Probability MAPE"])
        y_max <- max(results_sum_prediction_c[,"Occupancy Probability MAPE"])


        slope_diag <- (y_max - y_min) / (x_max - x_min)
        int_diag <- y_max - (slope_diag * x_max)

        # Correlation test
        correlation_result <- cor.test(results_sum_prediction_c[,metrics_list_formatted[[metric_name]]], results_sum_prediction_c[,"Occupancy Probability MAPE"], method = "pearson") 
        pearson_r <- format(round(correlation_result$estimate, 2), nsmall = 2)
        p_value <- format(round(correlation_result$p.value, 2), nsmall = 2)

        annotate_hjust <- -2.7
        if (pearson_r<0){
            annotate_hjust <- -2.45
        }

        pl <- ggplot(results_sum_prediction_c, aes(x=!!as.name(metrics_list_formatted[[metric_name]]), y=!!as.name("Occupancy Probability MAPE"), color = method)) +
            geom_abline(
                slope = slope_diag,
                intercept = int_diag,
                color="#A9A9A9",
                linetype = "dashed"
            ) +
            geom_point(size=3, stroke=1, shape = 8) +
            theme_linedraw() +
            scale_color_manual(name = "Algorithm", labels = alg_o, values=colors_c, drop = FALSE, na.translate = FALSE) +
            labs(title = spec, subtitle = species_descr_ext[species_descr_ext$species==spec,]$species_ext) +
            theme(
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(), 
                plot.title = element_text(size = 11, hjust=0.5, face = "bold"),
                plot.subtitle = element_text(size = 10, hjust = 0.5)
            ) +
            annotate(geom="text", x=-Inf, y=-Inf, label=paste0("r = ", pearson_r), hjust = annotate_hjust, vjust = -2, color="black") +
            annotate(geom="text", x=-Inf, y=-Inf, label=paste0("p = ", p_value), hjust = -2.5, vjust = -0.45, color="black")
            
            
        
        
        plots[[i]] <- pl

        i <- i + 1

    }
    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=ncol_c, common.legend=TRUE, legend="right")

    plot_main <- annotate_figure(
        plot_main,
        left = text_grob("Negative of Occupancy Probability MAPE (↑)", size = 13, hjust = 0.5, color = "black", rot = 90),
        bottom = text_grob(sup_x_title, size = 13,  hjust = sup_x_title_hjust, ,color = "black")
    )    


    show(plot_main)
    dev.off()




}


############################################################################################
# 4. Cluster diameter, and area historgrams
############################################################################################



# rep_species <- "Sim-1"

# # Create a list to store the metrics for each method
# all_site_metrics <- list()


# for (m in alg_all_o_cl_hist){
#     # Construct file path and read the clustering results
#     # In a real run, you would create dummy files for this example to work.
#     # For now, we will simulate the pts_df data frame.
#     file_path <- paste0("results/simulation_experiments/raw/clusterings/", rep_species,"_", m, ".csv")

#     # Check if the file exists before trying to read it
#     pts_df <- read.csv(file_path, header = T)

#     # Calculate metrics for each site using dplyr and geosphere
#     site_metrics <- pts_df %>%
#         mutate(site = as.character(site)) %>% # Ensure site ID is character type for consistency
#         group_by(site) %>%
#         reframe(
#             n_points = n(),
#             diameter_m = {
#                 coords <- cbind(longitude, latitude)
#                 if (n() < 2) {
#                     0
#                 } else {
#                     max(distm(coords, fun = distHaversine))
#                 }
#             },
#             area_sq_m = {
#                 coords <- cbind(longitude, latitude)
#                 if (n() < 3) {
#                     0
#                 } else {
#                     tryCatch({
#                         hull_indices <- chull(coords)
#                         hull_coords <- coords[hull_indices, , drop = FALSE]
#                         areaPolygon(hull_coords)
#                     }, error = function(e) {
#                         return(0)
#                     })
#                 }
#             }
#         ) %>%
#         mutate(method = m)

#     # Store the results for this method in our list
#     all_site_metrics[[m]] <- site_metrics


# }

# # Combine the list of data frames into a single data frame
# final_metrics_df <- bind_rows(all_site_metrics)


# write.csv(final_metrics_df, "results/simulation_experiments/summarized/area_diameter_metrics.csv", row.names=FALSE)

# ############################################################################################
# # New section: Calculate summary statistics before plotting
# ############################################################################################

# # Calculate statistics for diameter (mean, median, min, max, and count of zero-diameter sites)
# diam_stats <- final_metrics_df %>%
#   group_by(method) %>%
#   summarize(
#     zero_count = sum(diameter_m == 0),
#     mean_val = mean(diameter_m[diameter_m > 0], na.rm = TRUE),
#     median_val = median(diameter_m[diameter_m > 0], na.rm = TRUE),
#     min_val = min(diameter_m[diameter_m > 0], na.rm = TRUE),
#     max_val = max(diameter_m[diameter_m > 0], na.rm = TRUE),
#     .groups = 'drop'
#   )

# # Calculate statistics for area (mean, median, min, max, and count of zero-area sites)
# area_stats <- final_metrics_df %>%
#   group_by(method) %>%
#   summarize(
#     zero_count = sum(area_sq_m == 0),
#     mean_val = mean(area_sq_m[area_sq_m > 0], na.rm = TRUE),
#     median_val = median(area_sq_m[area_sq_m > 0], na.rm = TRUE),
#     min_val = min(area_sq_m[area_sq_m > 0], na.rm = TRUE),
#     max_val = max(area_sq_m[area_sq_m > 0], na.rm = TRUE),
#     .groups = 'drop'
#   )

# ############################################################################################
# # Updated plotting section with pre-filtering
# ############################################################################################



# # Create a filtered dataset that excludes all zero-value sites
# # This will be the source for all histograms.
# metrics_positive_only <- final_metrics_df %>%
#   filter(diameter_m > 0 | area_sq_m > 0)

# create_hist_plot <- function(data, method_name, metric, stats_df) {
#   # Filter data for the specific method
#   plot_data <- data %>% filter(method == method_name)

#   # Get the pre-calculated stats for this method
#   stats <- stats_df %>% filter(method == method_name)

#   # Define axis label and clean name for the metric
#   x_label <- if (metric == "diameter_m") "Diameter (m)" else "Area (m²)"
#   unit <- if (metric == "diameter_m") "m" else "m²"
#   metric_clean_name <- if (metric == "diameter_m") "Diameter" else "Area"

#   y_label <- "Count"

#   # Define statistics and colors
#   stat_names <- c("Min", "Mean", "Median", "Max")
#   stat_values <- c(stats$min_val, stats$mean_val, stats$median_val, stats$max_val)
#   stat_colors <- c("darkgreen", "red", "blue", "darkorange")
  
#   # Create a data frame for the vlines
#   vline_data <- data.frame(
#     name = factor(stat_names, levels = stat_names),
#     value = stat_values,
#     color = stat_colors
#   )

#   # --- STEP 1: Create the base plot ---
#   p_base <- ggplot(plot_data, aes(x = .data[[metric]])) +
#     theme_bw() +
#     geom_histogram(bins = 30, fill = "lightyellow", color = "black")
  
#   # --- STEP 2: Calculate Y-positions for the labels ---
#   max_y <- ggplot_build(p_base)$layout$panel_params[[1]]$y.range[2]
  
#   # Create a new data frame for the labels with staggered y-positions
#   # REVISED LOGIC: Sort by the stat's value before assigning Y positions.
#   # This prevents label connector lines from crossing over each other.
#   label_data <- vline_data %>%
#     arrange(value) %>% # <-- NEW: Sort by the statistic's value
#     mutate(
#       label_text = paste0(name, ": ", round(value, 2)),
#       y_pos = max_y * seq(0.95, 0.70, length.out = n()) # Stagger from 95% down to 70%
#     )

#   # --- STEP 3: Add all layers to the final plot ---
#   p_final <- p_base +
#     # Add the vertical lines
#     geom_vline(
#         data = vline_data,
#         aes(xintercept = value, color = name),
#         linetype = "dashed",
#         linewidth = 1
#     ) +
#     # Add the text labels next to the lines using geom_label_repel
#     geom_label_repel(
#         data = label_data,
#         aes(x = value, y = y_pos, label = label_text, color = name),
#         direction = "y",      # Force labels to move vertically
#         nudge_x = 0.1 * max(plot_data[[metric]]), 
#         segment.size = 0.5,
#         # --- NEW PARAMETERS TO PREVENT OVERLAP ---
#         box.padding = 0.5,         # Increase padding around the label box
#         force = 30,                # Increase the repulsion force
#         min.segment.length = 0,    # Always draw segment
#         max.overlaps = Inf,        # Ensure all labels are shown
#         seed = 42
#     ) +
#     # Set the colors for both lines and labels
#     scale_color_manual(
#         values = setNames(stat_colors, stat_names)
#     ) +
#     # Add titles and labels
#     ggtitle(paste0(method_name, "\n(Sites with 0 ", unit, " ", tolower(metric_clean_name), ": ", stats$zero_count, ")")) +
#     xlab(x_label) +
#     ylab(y_label) +
#     # --- STEP 4: Hide the legend ---
#     theme(legend.position = "none")

#   return(p_final)
# }

# # Set plot dimensions
# width_p <- 10
# height_p <- 8

# # Plot site diameter histograms using the pre-filtered data
# png("results/simulation_experiments/plotted/plots/diameter_hist.png", units="in", width=width_p, height=height_p, res=300)
# p1 <- create_hist_plot(metrics_positive_only, "1-kmSq", "diameter_m", diam_stats)
# p2 <- create_hist_plot(metrics_positive_only, "rounded-4", "diameter_m", diam_stats)
# p3 <- create_hist_plot(metrics_positive_only, "DBSC", "diameter_m", diam_stats)
# p4 <- create_hist_plot(metrics_positive_only, "BayesOptClustGeo", "diameter_m", diam_stats)
# print((p1 | p2) / (p3 | p4) + plot_annotation(title = "Distribution of Site Diameters (sites > 0m)"))
# dev.off()



# # Set plot dimensions
# width_p <- 15
# height_p <- 4

# # Plot site area histograms using the pre-filtered data
# png("results/simulation_experiments/plotted/plots/area_hist.png", units="in", width=width_p, height=height_p, res=300)
# p1 <- create_hist_plot(metrics_positive_only, "1-kmSq", "area_sq_m", area_stats)
# p2 <- create_hist_plot(metrics_positive_only, "DBSC", "area_sq_m", area_stats)
# p3 <- create_hist_plot(metrics_positive_only, "BayesOptClustGeo", "area_sq_m", area_stats)
# print((p1 | p2 | p3) + plot_annotation(title = "Distribution of Site Areas (sites > 0m²)"))
# dev.off()


rep_species <- "Sim-1"

# Create a list to store the metrics for each method
all_site_metrics <- list()


for (m in alg_all_o_cl_hist){
    # Construct file path and read the clustering results
    # In a real run, you would create dummy files for this example to work.
    # For now, we will simulate the pts_df data frame.
    file_path <- paste0("results/simulation_experiments/raw/clusterings/", rep_species,"_", m, ".csv")

    # Check if the file exists before trying to read it
    pts_df <- read.csv(file_path, header = T)

    # Calculate metrics for each site using dplyr and geosphere
    site_metrics <- pts_df %>%
        mutate(site = as.character(site)) %>% # Ensure site ID is character type for consistency
        group_by(site) %>%
        reframe(
            n_points = n(),
            diameter_m = {
                coords <- cbind(longitude, latitude)
                if (n() < 2) {
                    0
                } else {
                    max(distm(coords, fun = distHaversine))
                }
            },
            area_sq_m = {
                coords <- cbind(longitude, latitude)
                if (n() < 3) {
                    0
                } else {
                    tryCatch({
                        hull_indices <- chull(coords)
                        hull_coords <- coords[hull_indices, , drop = FALSE]
                        areaPolygon(hull_coords)
                    }, error = function(e) {
                        return(0)
                    })
                }
            }
        ) %>%
        mutate(method = m)

    # Store the results for this method in our list
    all_site_metrics[[m]] <- site_metrics


}

# Combine the list of data frames into a single data frame
final_metrics_df <- bind_rows(all_site_metrics)


write.csv(final_metrics_df, "results/simulation_experiments/summarized/area_diameter_metrics.csv", row.names=FALSE)

############################################################################################
# New section: Conditionally convert units and calculate summary statistics
############################################################################################

# FLAG: Set to TRUE to plot in kilometers, FALSE to plot in meters.
in_km <- TRUE

if (in_km) {
  final_metrics_df <- final_metrics_df %>%
    mutate(
      diameter_km = diameter_m / 1000,
      area_sq_km = area_sq_m / 1000000
    )
  diam_metric_col <- "diameter_km"
  area_metric_col <- "area_sq_km"
  diam_title_unit <- "km"
  area_title_unit <- "km²"
} else {
  diam_metric_col <- "diameter_m"
  area_metric_col <- "area_sq_m"
  diam_title_unit <- "m"
  area_title_unit <- "m²"
}


# Calculate statistics for diameter using the selected unit
diam_stats <- final_metrics_df %>%
  group_by(method) %>%
  summarize(
    zero_count = sum(.data[[diam_metric_col]] == 0),
    mean_val = mean(.data[[diam_metric_col]][.data[[diam_metric_col]] > 0], na.rm = TRUE),
    median_val = median(.data[[diam_metric_col]][.data[[diam_metric_col]] > 0], na.rm = TRUE),
    min_val = min(.data[[diam_metric_col]][.data[[diam_metric_col]] > 0], na.rm = TRUE),
    max_val = max(.data[[diam_metric_col]][.data[[diam_metric_col]] > 0], na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate statistics for area using the selected unit
area_stats <- final_metrics_df %>%
  group_by(method) %>%
  summarize(
    zero_count = sum(.data[[area_metric_col]] == 0),
    mean_val = mean(.data[[area_metric_col]][.data[[area_metric_col]] > 0], na.rm = TRUE),
    median_val = median(.data[[area_metric_col]][.data[[area_metric_col]] > 0], na.rm = TRUE),
    min_val = min(.data[[area_metric_col]][.data[[area_metric_col]] > 0], na.rm = TRUE),
    max_val = max(.data[[area_metric_col]][.data[[area_metric_col]] > 0], na.rm = TRUE),
    .groups = 'drop'
  )

############################################################################################
# Updated plotting section with dynamic units and pre-filtering
############################################################################################

# Create a filtered dataset that excludes all zero-value sites
metrics_positive_only <- final_metrics_df %>%
  filter(.data[[diam_metric_col]] > 0 | .data[[area_metric_col]] > 0)

create_hist_plot <- function(data, method_name, metric, stats_df) {
  # Filter data for the specific method
  plot_data <- data %>% filter(method == method_name)

  # Get the pre-calculated stats for this method
  stats <- stats_df %>% filter(method == method_name)
  
  # Dynamically determine labels and units based on the metric column name
  is_km <- grepl("_km", metric, fixed = TRUE)
  is_diameter <- grepl("diameter", metric, fixed = TRUE)

  if (is_diameter) {
      x_label <- if (is_km) "Diameter (km)" else "Diameter (m)"
      unit <- if (is_km) "km" else "m"
      metric_clean_name <- "Diameter"
  } else {
      x_label <- if (is_km) "Area (km²)" else "Area (m²)"
      unit <- if (is_km) "km²" else "m²"
      metric_clean_name <- "Area"
  }
  y_label <- "Count"

  # Define statistics and colors
  stat_names <- c("Min", "Mean", "Median", "Max")
  stat_values <- c(stats$min_val, stats$mean_val, stats$median_val, stats$max_val)
  stat_colors <- c("darkgreen", "red", "blue", "darkorange")
  
  # Create a data frame for the vlines
  vline_data <- data.frame(
    name = factor(stat_names, levels = stat_names),
    value = stat_values,
    color = stat_colors
  )

  # --- STEP 1: Create the base plot ---
  p_base <- ggplot(plot_data, aes(x = .data[[metric]])) +
    theme_bw() +
    geom_histogram(bins = 30, fill = "lightyellow", color = "black")
  
  # --- STEP 2: Calculate Y-positions for the labels ---
  max_y <- ggplot_build(p_base)$layout$panel_params[[1]]$y.range[2]
  
  # Create a new data frame for the labels with staggered y-positions
  label_data <- vline_data %>%
    arrange(value) %>% 
    mutate(
      label_text = sprintf("%s: %.3g", name, value), # UPDATED LINE
      y_pos = max_y * seq(0.95, 0.70, length.out = n()) # Stagger from 95% down to 70%
    )

  # --- STEP 3: Add all layers to the final plot ---
  p_final <- p_base +
    # Add the vertical lines
    geom_vline(
        data = vline_data,
        aes(xintercept = value, color = name),
        linetype = "dashed",
        linewidth = 1
    ) +
    # Add the text labels next to the lines using geom_label_repel
    geom_label_repel(
        data = label_data,
        aes(x = value, y = y_pos, label = label_text, color = name),
        direction = "y",
        nudge_x = 0.1 * max(plot_data[[metric]]), 
        segment.size = 0.5,
        box.padding = 0.5,
        force = 30,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 42
    ) +
    # Set the colors for both lines and labels
    scale_color_manual(
        values = setNames(stat_colors, stat_names)
    ) +
    # Add titles and labels
    ggtitle(paste0(method_name, "\n(Sites with 0 ", unit, " ", tolower(metric_clean_name), ": ", stats$zero_count, ")")) +
    xlab(x_label) +
    ylab(y_label) +
    # --- STEP 4: Hide the legend ---
    theme(legend.position = "none")

  return(p_final)
}

# Set plot dimensions
width_p <- 10
height_p <- 8

# Plot site diameter histograms using the selected unit
png(paste0("results/simulation_experiments/plotted/plots/diameter_hist_", diam_title_unit, ".png"), units="in", width=width_p, height=height_p, res=300)
p1 <- create_hist_plot(metrics_positive_only, "1-kmSq", diam_metric_col, diam_stats)
p2 <- create_hist_plot(metrics_positive_only, "rounded-4", diam_metric_col, diam_stats)
p3 <- create_hist_plot(metrics_positive_only, "DBSC", diam_metric_col, diam_stats)
p4 <- create_hist_plot(metrics_positive_only, "BayesOptClustGeo", diam_metric_col, diam_stats)
print((p1 | p2) / (p3 | p4) + plot_annotation(title = paste0("Distribution of Site Diameters (sites > 0", diam_title_unit, ")")))
dev.off()



# Set plot dimensions
width_p <- 15
height_p <- 4

# Plot site area histograms using the selected unit
png(paste0("results/simulation_experiments/plotted/plots/area_hist_", gsub("²", "2", area_title_unit), ".png"), units="in", width=width_p, height=height_p, res=300)
p1 <- create_hist_plot(metrics_positive_only, "1-kmSq", area_metric_col, area_stats)
p2 <- create_hist_plot(metrics_positive_only, "DBSC", area_metric_col, area_stats)
p3 <- create_hist_plot(metrics_positive_only, "BayesOptClustGeo", area_metric_col, area_stats)
print((p1 | p2 | p3) + plot_annotation(title = paste0("Distribution of Site Areas (sites > 0", area_title_unit, ")")))
dev.off()