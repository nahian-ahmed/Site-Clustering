library(ggplot2)
library(terra)
library(sf)
library(dplyr)
library(ggforce) # For geom_mark_ellipse
library(patchwork) # For combining plots
library(ggnewscale) # For multiple color scales


plot_sites <- function(
    base_train_df,
    all_clusterings,
    all_site_geometries, # <-- NEW ARGUMENT
    elevation_raster,
    methods_to_plot,
    boundary_shp_path,
    output_path,
    zoom_box = list(
        longitude = c(-123.025, -122.992),
        latitude = c(44.085, 44.118)
    )
) {

    # --- 1. Prepare Base Raster for Plotting ---
    valid_boundary <- terra::vect(boundary_shp_path)
    terra::crs(valid_boundary) <- terra::crs(elevation_raster)
    
    # Use the first layer (elevation)
    region <- terra::crop(elevation_raster[[1]], valid_boundary, mask = TRUE)
    
    base_rast_df <- as.data.frame(region, xy = TRUE)
    # Get the name of the elevation column
    elev_col_name <- names(base_rast_df)[3] 
    
    base_rast_min <- min(base_rast_df[[elev_col_name]], na.rm = TRUE)
    base_rast_max <- max(base_rast_df[[elev_col_name]], na.rm = TRUE)

    bbox_full <- terra::ext(region)
    
    # --- Define WGS84 CRS for transformations ---
    wgs84_crs <- "+proj=longlat +datum=WGS84"

    # --- 2. Create Left Plot (Observations) ---
    obs_plot <- ggplot() +
        geom_raster(
            data = base_rast_df, 
            aes(x = x, y = y, fill = .data[[elev_col_name]]), 
            show.legend = TRUE
        ) +
        scale_fill_viridis_c(
            option = "H", 
            limits = c(base_rast_min, base_rast_max),
            name = "Elevation (m)"
        ) +
        geom_point(
            data = base_train_df, 
            aes(x = longitude, y = latitude), 
            color = "black", size = 1, shape = 21, fill = "#FBFAF5", 
            show.legend = FALSE
        ) +
        annotate(
            "segment",
            x = -123.75, xend = zoom_box$longitude[1],
            y = 44.425, yend = (zoom_box$latitude[1] + zoom_box$latitude[2]) / 2,
            colour = "yellow"
        ) +
        annotate(
            "label",
            x = -124.05, y = 44.425,
            label = "Zoomed-in region",
            size = 3.25,
        ) +
        geom_rect(
            aes(
                xmin = zoom_box$longitude[1], xmax = zoom_box$longitude[2],
                ymin = zoom_box$latitude[1], ymax = zoom_box$latitude[2]
            ), 
            fill = "red", color = "red", alpha = 0.1, show.legend = FALSE
        ) +
        labs(
            title = "Species Observations",
            x = "Longitude",
            y = "Latitude"
        ) +
        theme_bw() +
        # Use coord_sf to set map limits and CRS
        coord_sf(
            xlim = c(bbox_full$xmin, bbox_full$xmax),
            ylim = c(bbox_full$ymin, bbox_full$ymax),
            crs = wgs84_crs,
            expand = FALSE
        ) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.width = unit(1.5, "cm"),
            axis.title = element_blank() # Remove axis titles for maps
        )

    # --- 3. Create Right Plots (Zoomed Clusters) ---
    
    zoom_plots <- list()
    
    for (method_name in methods_to_plot) {
        
        # --- Get Point Data ---
        if (!method_name %in% names(all_clusterings)) {
            warning(paste("Method not found in all_clusterings:", method_name))
            next
        }
        pts_df <- all_clusterings[[method_name]]
        if (is.list(pts_df) && "result_df" %in% names(pts_df)) {
          pts_df <- pts_df$result_df
        }
        
        # --- Get Geometry Data ---
        if (!method_name %in% names(all_site_geometries)) {
            warning(paste("Method not found in all_site_geometries:", method_name))
            next
        }
        geom_sf <- all_site_geometries[[method_name]]
        if (is.null(geom_sf)) {
            warning(paste("NULL geometry for method:", method_name))
            next
        }
        
        # --- Transform Geometries to WGS84 for plotting ---
        geom_sf_wgs84 <- sf::st_transform(geom_sf, crs = wgs84_crs)
        geom_sf_wgs84$site <- as.factor(geom_sf_wgs84$site)
        

        # --- Filter Point Data ---
        pts_df$site <- as.factor(pts_df$site)
        pts_df_zoom <- pts_df[
            (pts_df$longitude > zoom_box$longitude[1]) &
            (pts_df$longitude < zoom_box$longitude[2]) &
            (pts_df$latitude > zoom_box$latitude[1]) &
            (pts_df$latitude < zoom_box$latitude[2]), ]
        
        pts_df_distinct <- pts_df_zoom %>%
            distinct(latitude, longitude, .keep_all = TRUE)
            
        p_zoom <- ggplot() +
            # Faded background
            geom_raster(
                data = base_rast_df, 
                aes(x = x, y = y, fill = .data[[elev_col_name]]),
                show.legend = FALSE
            ) +
            scale_fill_viridis_c(
                option = "H",
                limits = c(base_rast_min, base_rast_max),
                aesthetics = "fill"
            ) +
            new_scale_fill() + # Allows for a new fill scale
            
            # --- NEW LAYER: Site Geometries ---
            geom_sf(
                data = geom_sf_wgs84,
                aes(fill = site), # Color by site
                alpha = 0.4,      # Make semi-transparent
                color = "black",  # Add a black border
                linewidth = 0.5,
                show.legend = FALSE
            ) +
            
            # --- Clustered points (ON TOP) ---
            geom_point(
                data = pts_df_distinct,
                aes(x = longitude, y = latitude, fill = site),
                shape = 21, size = 2.5, # Slightly smaller points
                color = "black",
                show.legend = FALSE
            ) +
            
            # --- REMOVED geom_mark_ellipse ---
            
            scale_fill_discrete() + # This scale applies to both geom_sf and geom_point
            
            # --- Use coord_sf to set zoom and CRS ---
            coord_sf(
                xlim = c(zoom_box$longitude[1], zoom_box$longitude[2]),
                ylim = c(zoom_box$latitude[1], zoom_box$latitude[2]),
                crs = wgs84_crs,
                expand = FALSE # Don't expand beyond the zoom box
            ) +
            theme_bw() +
            labs(title = method_name) +
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 10, hjust = 0.5),
                panel.grid.major = element_blank(), # Remove gridlines
                panel.grid.minor = element_blank()
            )
            
        zoom_plots[[method_name]] <- p_zoom
    }
    
    # --- 4. Assemble the Final Plot ---
    
    # Arrange the grid of zoom plots
    plot_clust <- patchwork::wrap_plots(zoom_plots, ncol = 3)
    
    # Combine the main plot and the grid
    final_plot <- obs_plot + plot_clust +
        plot_layout(nrow = 1, widths = c(1, 1.5)) # Adjust width ratio as needed



    ggsave(
        output_path,
        plot = final_plot,
        width = 16,
        height = 8,
        dpi = 300
    )
    cat(sprintf("--- Site cluster plot saved to %s/site_cluster_visualization.png ---\n", output_dir))

    
}