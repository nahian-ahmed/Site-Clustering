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
    ),
    cluster_labels = FALSE
) {

    # --- 0. Setup Coordinate Systems ---
    albers_crs_str <- terra::crs(elevation_raster) 
    wgs84_crs_sf <- sf::st_crs(4326)               
    wgs84_crs_str <- "EPSG:4326"                   

    # --- 1. Prepare Rasters ---
    
    # A. Prepare ALBERS Raster
    valid_boundary <- terra::vect(boundary_shp_path)
    valid_boundary_albers <- terra::project(valid_boundary, albers_crs_str)
    
    region_albers <- terra::crop(elevation_raster[[1]], valid_boundary_albers, mask = TRUE)
    base_rast_df_albers <- as.data.frame(region_albers, xy = TRUE)
    elev_col_name <- names(base_rast_df_albers)[3]
    
    base_rast_min <- min(base_rast_df_albers[[elev_col_name]], na.rm = TRUE)
    base_rast_max <- max(base_rast_df_albers[[elev_col_name]], na.rm = TRUE)

    # B. Prepare WGS84 Raster
    region_wgs84 <- terra::project(region_albers, wgs84_crs_str)
    base_rast_df_wgs84 <- as.data.frame(region_wgs84, xy = TRUE)
    bbox_full <- terra::ext(region_wgs84) 

    
    # --- 2. Create Left Plot (Observations in WGS84) ---
    obs_plot <- ggplot() +
        geom_raster(
            data = base_rast_df_wgs84,
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
            x = -123.85, xend = zoom_box$longitude[1],
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
        coord_fixed(
            ratio = 1.0,
            xlim = c(bbox_full$xmin, bbox_full$xmax),
            ylim = c(bbox_full$ymin, bbox_full$ymax),
            expand = FALSE
        ) +
        theme(
            plot.title = element_text(
                size = 12,
                hjust = 0.5,
                face = "bold",
                vjust = -1.5, 
                margin = margin(b = -10)
            ),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.margin = margin(t = -50),
            legend.box.margin = margin(0, 0, 0, 0),
            legend.key.width = unit(0.75, "cm"),
            legend.key.height = unit(0.5, "cm"),
            legend.title = element_text(size = 10, vjust = 1),
            legend.text = element_text(size = 8),
            axis.title = element_text(size = 10)
        )

    # --- 3. Create Right Plots (Zoomed Clusters in ALBERS) ---
    
    zoom_poly_wgs84 <- sf::st_as_sfc(sf::st_bbox(c(
        xmin = zoom_box$longitude[1], xmax = zoom_box$longitude[2],
        ymin = zoom_box$latitude[1], ymax = zoom_box$latitude[2]
    ), crs = wgs84_crs_sf)) 
    
    zoom_poly_albers <- sf::st_transform(zoom_poly_wgs84, albers_crs_str)
    zoom_bbox_albers <- sf::st_bbox(zoom_poly_albers) 

    zoom_plots <- list()

    for (method_name in methods_to_plot) {

        # --- Get Point Data & Transform to Albers ---
        if (!method_name %in% names(all_clusterings)) next
        pts_df <- all_clusterings[[method_name]]
        if (is.list(pts_df) && "result_df" %in% names(pts_df)) pts_df <- pts_df$result_df
        
        # Project points to Albers
        pts_sf <- sf::st_as_sf(pts_df, coords = c("longitude", "latitude"), crs = wgs84_crs_sf)
        pts_sf_albers <- sf::st_transform(pts_sf, albers_crs_str)
        
        # Extract coordinates for plotting
        pts_coords <- sf::st_coordinates(pts_sf_albers)
        pts_df_albers <- pts_df
        pts_df_albers$x <- pts_coords[,1]
        pts_df_albers$y <- pts_coords[,2]
        pts_df_albers$site <- as.factor(pts_df_albers$site)

        # Filter points to zoom box
        pts_df_zoom <- pts_df_albers[
            (pts_df_albers$x > zoom_bbox_albers["xmin"]) &
            (pts_df_albers$x < zoom_bbox_albers["xmax"]) &
            (pts_df_albers$y > zoom_bbox_albers["ymin"]) &
            (pts_df_albers$y < zoom_bbox_albers["ymax"]), 
        ]
        
        # --- Pre-calculate Labels (1..N) and Position if requested ---
        if (cluster_labels && nrow(pts_df_zoom) > 0) {
            # 1. Map Site IDs to sequential 1..N based on factor order
            # This ensures sites are labeled 1, 2, 3... within this specific subplot
            visible_sites <- levels(droplevels(pts_df_zoom$site))
            site_map <- seq_along(visible_sites)
            names(site_map) <- visible_sites
            
            pts_df_zoom$plot_label <- site_map[as.character(pts_df_zoom$site)]
            
            # 2. Determine HJUST based on proximity to right margin
            # If x is in the top 15% of the range, flip label to left
            x_range <- zoom_bbox_albers["xmax"] - zoom_bbox_albers["xmin"]
            right_margin_thresh <- zoom_bbox_albers["xmax"] - (0.15 * x_range)
            
            # -0.4 places text to the right (default), 1.2 places text to the left (if near margin)
            pts_df_zoom$lab_hjust <- ifelse(pts_df_zoom$x > right_margin_thresh, 1.4, -0.4)
        }

        # --- Get Geometry Data (Already Albers) ---
        if (!method_name %in% names(all_site_geometries)) next
        geom_sf <- all_site_geometries[[method_name]]
        if (is.null(geom_sf)) next
        
        if (is.na(sf::st_crs(geom_sf))) sf::st_crs(geom_sf) <- albers_crs_str
        geom_sf$site <- as.factor(geom_sf$site)

        # Crop geometries
        geom_sf_zoom <- suppressWarnings(sf::st_crop(geom_sf, zoom_bbox_albers))
        
        if (nrow(geom_sf_zoom) > 0) {
             if (any(grepl("COLLECTION", sf::st_geometry_type(geom_sf_zoom)))) {
                 geom_sf_zoom <- sf::st_collection_extract(geom_sf_zoom, "POLYGON")
             }
        }

        # --- Plotting (Albers) ---
        p_zoom <- ggplot() +
            # Use Albers raster dataframe
            geom_raster(
                data = base_rast_df_albers,
                aes(x = x, y = y, fill = .data[[elev_col_name]]),
                show.legend = FALSE
            ) +
            scale_fill_viridis_c(
                option = "H",
                limits = c(base_rast_min, base_rast_max),
                aesthetics = "fill"
            ) +
            new_scale_fill() +

            # Site Geometries (Albers)
            {if (nrow(geom_sf_zoom) > 0)
                geom_sf(
                    data = geom_sf_zoom,
                    aes(fill = site),
                    alpha = 0.4,
                    color = "black",
                    linewidth = 0.25,
                    show.legend = FALSE,
                    inherit.aes = FALSE
                )
            } +

            # Clustered points (Albers)
            geom_point(
                data = pts_df_zoom,
                aes(x = x, y = y, fill = site),
                shape = 21, size = 2.0,
                color = "black",
                show.legend = FALSE
            ) +
            
            # --- NEW: Conditional Labels ---
            {if (cluster_labels && nrow(pts_df_zoom) > 0)
                geom_text(
                    data = pts_df_zoom,
                    aes(
                        x = x, 
                        y = y, 
                        label = plot_label, 
                        hjust = lab_hjust
                    ),
                    vjust = -0.5, # Move up slightly (exponential style)
                    size = 2.5,   # Small text
                    fontface = "bold",
                    color = "black",
                    inherit.aes = FALSE
                )
            } +
            
            scale_fill_discrete() +
            
            # Enforce Albers Coordinates
            coord_sf(
                xlim = c(zoom_bbox_albers["xmin"], zoom_bbox_albers["xmax"]),
                ylim = c(zoom_bbox_albers["ymin"], zoom_bbox_albers["ymax"]),
                crs = albers_crs_str,
                expand = FALSE
            ) +
            theme_bw() +
            labs(title = method_name) +
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 10, hjust = 0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
            )

        zoom_plots[[method_name]] <- p_zoom
    }

    # --- 4. Assemble the Final Plot ---
    plot_clust <- patchwork::wrap_plots(zoom_plots, ncol = 6)
    final_plot <- obs_plot + plot_clust +
        plot_layout(nrow = 1, widths = c(1, 4))

    ggsave(
        output_path,
        plot = final_plot,
        width = 14,
        height = 6, 
        dpi = 300
    )

    cat(sprintf("--- Site cluster plot saved to %s ---\n", output_path))
}