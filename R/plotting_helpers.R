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
    elevation_raster, # Assumed to be the Albers Raster
    methods_to_plot,
    boundary_shp_path,
    output_path,
    zoom_box = list(
        longitude = c(-123.025, -122.992),
        latitude = c(44.085, 44.118)
    )
) {

    cat("--- Generating Site Plots (Projecting to Raster CRS) ---\n")

    # 1. Detect Raster CRS (Albers)
    target_crs <- terra::crs(elevation_raster)
    
    # 2. Prepare Base Raster DataFrame
    valid_boundary <- terra::vect(boundary_shp_path)
    valid_boundary <- terra::project(valid_boundary, target_crs)
    
    # Crop to boundary to reduce plot size
    region <- terra::crop(elevation_raster[[1]], valid_boundary, mask = TRUE)
    base_rast_df <- as.data.frame(region, xy = TRUE)
    elev_col_name <- names(base_rast_df)[3] 
    
    # 3. Project Training Points to Match Raster
    points_sf <- sf::st_as_sf(
        base_train_df, 
        coords = c("longitude", "latitude"), 
        crs = 4326 # WGS84
    )
    points_proj <- sf::st_transform(points_sf, target_crs)
    
    # --- FIX: Safely add coordinates without duplicate column names ---
    points_df_proj <- sf::st_drop_geometry(points_proj)
    coords <- sf::st_coordinates(points_proj)
    points_df_proj$x <- coords[, 1]
    points_df_proj$y <- coords[, 2]
    # ----------------------------------------------------------------

    # 4. Project Zoom Box
    zoom_poly_wgs84 <- sf::st_polygon(list(rbind(
        c(zoom_box$longitude[1], zoom_box$latitude[1]),
        c(zoom_box$longitude[2], zoom_box$latitude[1]),
        c(zoom_box$longitude[2], zoom_box$latitude[2]),
        c(zoom_box$longitude[1], zoom_box$latitude[2]),
        c(zoom_box$longitude[1], zoom_box$latitude[1])
    )))
    zoom_sfc_wgs84 <- sf::st_sfc(zoom_poly_wgs84, crs = 4326)
    zoom_sfc_proj <- sf::st_transform(zoom_sfc_wgs84, target_crs)
    zoom_bbox <- sf::st_bbox(zoom_sfc_proj) 

    # --- LEFT PLOT: OBSERVATIONS ---
    obs_plot <- ggplot() +
        geom_raster(
            data = base_rast_df, 
            aes(x = x, y = y, fill = .data[[elev_col_name]])
        ) +
        scale_fill_viridis_c(
            option = "H", 
            name = "Elevation (m)"
        ) +
        geom_point(
            data = points_df_proj, 
            aes(x = x, y = y), 
            color = "black", size = 0.8, shape = 21, fill = "#FBFAF5"
        ) +
        geom_sf(data = zoom_sfc_proj, fill = NA, color = "red", linewidth = 0.8) +
        annotate(
            "label",
            x = zoom_bbox["xmin"], y = zoom_bbox["ymax"] + 200, 
            label = "Zoomed Region",
            size = 3, fontface = "bold"
        ) +
        coord_sf(expand = FALSE) + 
        theme_bw() +
        labs(title = "Species Observations (Albers Projection)", x = "x (m)", y = "y (m)") +
        theme(
            legend.position = "bottom",
            legend.key.width = unit(1.5, "cm")
        )

    # --- RIGHT PLOTS: ZOOMED CLUSTERS ---
    zoom_plots <- list()
    
    for (method_name in methods_to_plot) {
        
        if (!method_name %in% names(all_clusterings)) next
        pts_df <- all_clusterings[[method_name]]
        if (is.list(pts_df) && "result_df" %in% names(pts_df)) pts_df <- pts_df$result_df
        
        if (!method_name %in% names(all_site_geometries)) next
        geom_sf <- all_site_geometries[[method_name]]
        if (is.null(geom_sf)) next
        
        if (sf::st_crs(geom_sf) != sf::st_crs(target_crs)) {
            geom_sf <- sf::st_transform(geom_sf, target_crs)
        }
        
        pts_sf_method <- sf::st_as_sf(pts_df, coords = c("longitude", "latitude"), crs = 4326)
        pts_sf_method <- sf::st_transform(pts_sf_method, target_crs)
        
        geom_sf_zoom <- suppressWarnings(sf::st_crop(geom_sf, zoom_bbox))
        pts_sf_zoom  <- suppressWarnings(sf::st_crop(pts_sf_method, zoom_bbox))
        
        if (nrow(geom_sf_zoom) > 0) {
             if (any(grepl("COLLECTION", sf::st_geometry_type(geom_sf_zoom)))) {
                 geom_sf_zoom <- sf::st_collection_extract(geom_sf_zoom, "POLYGON")
             }
        }

        p_zoom <- ggplot() +
            geom_raster(
                data = base_rast_df, 
                aes(x = x, y = y, fill = .data[[elev_col_name]]),
                show.legend = FALSE
            ) +
            scale_fill_viridis_c(option = "H") +
            new_scale_fill() + 
            geom_sf(
                data = geom_sf_zoom, 
                aes(fill = as.factor(site)), 
                alpha = 0.4, color = "black", linewidth = 0.2, show.legend = FALSE
            ) +
            geom_sf(
                data = pts_sf_zoom,
                aes(fill = as.factor(site)),
                shape = 21, size = 2.0, color = "black", show.legend = FALSE
            ) +
            coord_sf(
                xlim = c(zoom_bbox["xmin"], zoom_bbox["xmax"]),
                ylim = c(zoom_bbox["ymin"], zoom_bbox["ymax"]),
                expand = FALSE
            ) +
            theme_bw() +
            labs(title = method_name) +
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(size = 9, hjust = 0.5)
            )
            
        zoom_plots[[method_name]] <- p_zoom
    }
    
    cat("--- Assembling final plot... ---\n")
    plot_clust <- patchwork::wrap_plots(zoom_plots, ncol = 6)
    
    final_plot <- obs_plot + plot_clust +
        plot_layout(nrow = 1, widths = c(1, 4)) 

    ggsave(output_path, plot = final_plot, width = 16, height = 8, dpi = 300)
    cat(sprintf("--- Saved to %s ---\n", output_path))
}