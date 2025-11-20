library(ggplot2)
library(terra)
library(sf)
library(dplyr)
library(ggforce) 
library(patchwork) 
library(ggnewscale) 

plot_sites <- function(
    base_train_df,
    all_clusterings,
    all_site_geometries, 
    elevation_raster,
    methods_to_plot,
    boundary_shp_path,
    output_path,
    zoom_box = list(
        longitude = c(-123.025, -122.992),
        latitude = c(44.085, 44.118)
    )
) {

    # --- 0. Aspect Ratio Correction (Fixing the "Squeeze") ---
    # Calculate the center of the box
    mid_lon <- mean(zoom_box$longitude)
    mid_lat <- mean(zoom_box$latitude)
    
    # Determine height in degrees
    lat_height <- zoom_box$latitude[2] - zoom_box$latitude[1]
    
    # Calculate required width in degrees to make it visually square
    # Adjustment factor: 1 / cos(latitude in radians)
    lon_width_correction <- 1 / cos(mid_lat * pi / 180)
    required_lon_width <- lat_height * lon_width_correction
    
    # Update the zoom_box longitude to be physically square
    zoom_box$longitude <- c(mid_lon - required_lon_width/2, mid_lon + required_lon_width/2)

    # --- 1. Prepare Base Raster for Plotting ---
    valid_boundary <- terra::vect(boundary_shp_path)
    terra::crs(valid_boundary) <- terra::crs(elevation_raster)
    
    # Use the first layer (elevation)
    region <- terra::crop(elevation_raster[[1]], valid_boundary, mask = TRUE)
    
    base_rast_df <- as.data.frame(region, xy = TRUE)
    elev_col_name <- names(base_rast_df)[3] 
    
    base_rast_min <- min(base_rast_df[[elev_col_name]], na.rm = TRUE)
    base_rast_max <- max(base_rast_df[[elev_col_name]], na.rm = TRUE)

    bbox_full <- terra::ext(region)
    
    # --- Define WGS84 CRS (EPSG:4326) ---
    wgs84_crs <- sf::st_crs(4326)

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
            axis.title = element_blank() 
        )

    # --- 3. Create Right Plots (Zoomed Clusters) ---
    
    zoom_plots <- list()
    
    for (method_name in methods_to_plot) {
        
        # --- Get Point Data ---
        if (!method_name %in% names(all_clusterings)) {
            next
        }
        pts_df <- all_clusterings[[method_name]]
        if (is.list(pts_df) && "result_df" %in% names(pts_df)) {
          pts_df <- pts_df$result_df
        }
        
        # --- Get Geometry Data ---
        if (!method_name %in% names(all_site_geometries)) {
            next
        }
        geom_sf <- all_site_geometries[[method_name]]
        if (is.null(geom_sf)) {
            next
        }
        
        # --- Transform Geometries to WGS84 for plotting ---
        if (is.na(sf::st_crs(geom_sf))) {
             sf::st_crs(geom_sf) <- 5070 
        }
        geom_sf_wgs84 <- sf::st_transform(geom_sf, crs = wgs84_crs)
        geom_sf_wgs84$site <- as.factor(geom_sf_wgs84$site)
        
        
        # --- Robust Intersection Logic ---
        
        # 1. Create bbox
        zoom_bbox_sf <- sf::st_bbox(c(
            xmin = zoom_box$longitude[1], 
            xmax = zoom_box$longitude[2],
            ymin = zoom_box$latitude[1], 
            ymax = zoom_box$latitude[2]
        ), crs = wgs84_crs)
        zoom_poly_sfc <- sf::st_as_sfc(zoom_bbox_sf)

        # 2. Make valid BEFORE intersection
        sf::st_agr(geom_sf_wgs84) = "constant"
        geom_sf_wgs84 <- sf::st_make_valid(geom_sf_wgs84)

        # 3. Intersect
        geom_sf_zoom <- suppressWarnings(
            sf::st_intersection(geom_sf_wgs84, zoom_poly_sfc)
        )
        
        # 4. Post-Intersection Extraction
        if (nrow(geom_sf_zoom) > 0) {
            geom_sf_zoom <- sf::st_collection_extract(geom_sf_zoom, "POLYGON")
            geom_sf_zoom <- sf::st_make_valid(geom_sf_zoom)
            
            if (nrow(geom_sf_zoom) > 0) {
                area_vals <- sf::st_area(geom_sf_zoom)
                geom_sf_zoom <- geom_sf_zoom[as.numeric(area_vals) > 1e-6, ]
            }
        }

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
            new_scale_fill() + 
            
            # --- Site Geometries ---
            {if (nrow(geom_sf_zoom) > 0) 
                geom_sf(
                    data = geom_sf_zoom, 
                    aes(fill = site), 
                    alpha = 0.4,      
                    color = "black",  
                    linewidth = 0.3, 
                    show.legend = FALSE,
                    inherit.aes = FALSE 
                )
            } +
            
            # --- Clustered points ---
            geom_point(
                data = pts_df_distinct,
                aes(x = longitude, y = latitude, fill = site),
                shape = 21, size = 1.0, 
                color = "black",
                show.legend = FALSE
            ) +
            
            scale_fill_discrete() + 
            
            coord_sf(
                xlim = c(zoom_box$longitude[1], zoom_box$longitude[2]),
                ylim = c(zoom_box$latitude[1], zoom_box$latitude[2]),
                crs = wgs84_crs,
                expand = FALSE 
            ) +
            theme_bw() +
            labs(title = method_name) +
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 9, hjust = 0.5), # Reduced title size
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.margin = margin(0, 0, 2, 0) # Removed margins (Top, Right, Bottom, Left)
            )
            
        zoom_plots[[method_name]] <- p_zoom
    }
    
    # --- 4. Assemble the Final Plot ---
    # Reduced grid to 5 columns to allow more width per plot
    plot_clust <- patchwork::wrap_plots(zoom_plots, ncol = 6) 
    
    final_plot <- obs_plot + plot_clust +
        plot_layout(nrow = 1, widths = c(1.5, 6)) 

    ggsave(
        output_path,
        plot = final_plot,
        width = 18, # Slightly narrower to force height/width balance
        height = 9,
        dpi = 300
    )
    
    cat(sprintf("--- Site cluster plot saved to %s ---\n", output_path))
}