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
    # UPDATED: Using coord_fixed(ratio = 1) to force square aspect ratio
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
        # REPLACED coord_sf with coord_fixed for 1:1 lat/long ratio
        coord_fixed(
            ratio = 1,
            xlim = c(bbox_full$xmin, bbox_full$xmax),
            ylim = c(bbox_full$ymin, bbox_full$ymax),
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
        if (is.na(sf::st_crs(geom_sf))) {
             sf::st_crs(geom_sf) <- 5070 
        }
        geom_sf_wgs84 <- sf::st_transform(geom_sf, crs = wgs84_crs)
        geom_sf_wgs84$site <- as.factor(geom_sf_wgs84$site)
        
        
        # +++ FIX START: Robust Intersection Logic +++
        
        # 1. Create bbox
        zoom_bbox_sf <- sf::st_bbox(c(
            xmin = zoom_box$longitude[1], 
            xmax = zoom_box$longitude[2],
            ymin = zoom_box$latitude[1], 
            ymax = zoom_box$latitude[2]
        ), crs = wgs84_crs)
        zoom_poly_sfc <- sf::st_as_sfc(zoom_bbox_sf)

        # 2. Make valid BEFORE intersection (Crucial for fine grids)
        sf::st_agr(geom_sf_wgs84) = "constant"
        geom_sf_wgs84 <- sf::st_make_valid(geom_sf_wgs84)

        # 3. Intersect
        geom_sf_zoom <- suppressWarnings(
            sf::st_intersection(geom_sf_wgs84, zoom_poly_sfc)
        )
        
        # 4. Post-Intersection Extraction
        geom_poly_df <- data.frame() # Initialize empty for logic check later
        
        if (nrow(geom_sf_zoom) > 0) {
            
            # Ensure we extract POLYGONS from any GEOMETRYCOLLECTIONS
            geom_sf_zoom <- sf::st_collection_extract(geom_sf_zoom, "POLYGON")
            
            # Make valid again after extraction
            geom_sf_zoom <- sf::st_make_valid(geom_sf_zoom)
            
            # Filter by area
            if (nrow(geom_sf_zoom) > 0) {
                area_vals <- sf::st_area(geom_sf_zoom)
                geom_sf_zoom <- geom_sf_zoom[as.numeric(area_vals) > 1e-6, ]
            }
            
            # UPDATED: Convert SF to Dataframe for geom_polygon usage (allows coord_fixed)
            if (nrow(geom_sf_zoom) > 0) {
                # Cast to polygons to ensure consistency
                geom_poly_cast <- sf::st_cast(geom_sf_zoom, "POLYGON")
                # Extract coordinates
                geom_coords <- sf::st_coordinates(geom_poly_cast) %>% as.data.frame()
                # Map L1 (polygon index) back to Site ID
                geom_coords$site <- geom_poly_cast$site[geom_coords$L1]
                geom_poly_df <- geom_coords
            }
        }
        # +++ FIX END +++


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
            
            # --- Site Geometries (UPDATED: using geom_polygon) ---
            {if (nrow(geom_poly_df) > 0) 
                geom_polygon(
                    data = geom_poly_df, 
                    aes(x = X, y = Y, group = interaction(L1, L2), fill = site), 
                    alpha = 0.4,      
                    color = "black",  
                    linewidth = 0.5,
                    show.legend = FALSE
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
            
            # UPDATED: Replaced coord_sf with coord_fixed(ratio=1)
            coord_fixed(
                ratio = 1,
                xlim = c(zoom_box$longitude[1], zoom_box$longitude[2]),
                ylim = c(zoom_box$latitude[1], zoom_box$latitude[2]),
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
        plot_layout(nrow = 1, widths = c(2, 6)) 

    ggsave(
        output_path,
        plot = final_plot,
        width = 20,
        height = 8,
        dpi = 300
    )
    
    cat(sprintf("--- Site cluster plot saved to %s ---\n", output_path))
}