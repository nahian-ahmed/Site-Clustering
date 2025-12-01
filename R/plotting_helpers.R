# library(ggplot2)
# library(terra)
# library(sf)
# library(dplyr)
# library(ggforce) 
# library(patchwork) 
# library(ggnewscale) 

# plot_sites <- function(
#     base_train_df,
#     all_clusterings,
#     all_site_geometries, 
#     elevation_raster,
#     methods_to_plot,
#     boundary_shp_path,
#     output_path,
#     zoom_box = list(
#         longitude = c(-123.025, -122.992),
#         latitude = c(44.085, 44.118)
#     )
# ) {

#     # --- 1. Prepare Base Raster for Plotting ---
#     valid_boundary <- terra::vect(boundary_shp_path)
#     terra::crs(valid_boundary) <- terra::crs(elevation_raster)
    
#     # Use the first layer (elevation)
#     region <- terra::crop(elevation_raster[[1]], valid_boundary, mask = TRUE)
    
#     base_rast_df <- as.data.frame(region, xy = TRUE)
#     elev_col_name <- names(base_rast_df)[3] 
    
#     base_rast_min <- min(base_rast_df[[elev_col_name]], na.rm = TRUE)
#     base_rast_max <- max(base_rast_df[[elev_col_name]], na.rm = TRUE)

#     bbox_full <- terra::ext(region)
    
#     # --- Define WGS84 CRS (EPSG:4326) ---
#     wgs84_crs <- sf::st_crs(4326)

#     # --- 2. Create Left Plot (Observations) ---
#     obs_plot <- ggplot() +
#         geom_raster(
#             data = base_rast_df, 
#             aes(x = x, y = y, fill = .data[[elev_col_name]]), 
#             show.legend = TRUE
#         ) +
#         scale_fill_viridis_c(
#             option = "H", 
#             limits = c(base_rast_min, base_rast_max),
#             name = "Elevation (m)"
#         ) +
#         geom_point(
#             data = base_train_df, 
#             aes(x = longitude, y = latitude), 
#             color = "black", size = 1, shape = 21, fill = "#FBFAF5", 
#             show.legend = FALSE
#         ) +
#         annotate(
#             "segment",
#             x = -123.85, xend = zoom_box$longitude[1],
#             y = 44.425, yend = (zoom_box$latitude[1] + zoom_box$latitude[2]) / 2,
#             colour = "yellow"
#         ) +
#         annotate(
#             "label",
#             x = -124.05, y = 44.425,
#             label = "Zoomed-in region",
#             size = 3.25,
#         ) +
#         geom_rect(
#             aes(
#                 xmin = zoom_box$longitude[1], xmax = zoom_box$longitude[2],
#                 ymin = zoom_box$latitude[1], ymax = zoom_box$latitude[2]
#             ), 
#             fill = "red", color = "red", alpha = 0.1, show.legend = FALSE
#         ) +
#         labs(
#             title = "Species Observations",
#             x = "Longitude",
#             y = "Latitude"
#         ) +
#         theme_bw() +
#         coord_fixed(
#             ratio = 1.0,
#             xlim = c(bbox_full$xmin, bbox_full$xmax),
#             ylim = c(bbox_full$ymin, bbox_full$ymax),
#             expand = FALSE
#         ) +
#         theme(
            
#             # --- 1. Pull Title Closer ---
#             # Negative bottom margin (b = -15) pulls the title down towards the map
#             # plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = -100)),
#             plot.title = element_text(
#                 hjust = 0.5, 
#                 face = "bold", 
#                 vjust = -5,  # (Tweak number until it lands right) # -5 for 3 row, -1.25 for 2 row
#                 margin = margin(b = -10) # Keep a small margin adjustment if needed
#             ),

#             # --- 2. Legend Position & Margins ---
#             legend.position = "bottom",
#             legend.direction = "horizontal",
#             # Negative top margin (t = -15) pulls the legend up towards the map
#             legend.margin = margin(t = -250), # -50 for 2 row, -250 for 3 row
#             # Remove extra box spacing
#             legend.box.margin = margin(0, 0, 0, 0),

#             # --- 3. Shrink Legend Keys & Text ---
#             # Make the color bar much shorter and thinner
#             legend.key.width = unit(0.75, "cm"),  
#             legend.key.height = unit(0.5, "cm"),
            
#             # Reduce text sizes
#             legend.title = element_text(size = 10, vjust = 1),
#             legend.text = element_text(size = 8),

#             axis.title = element_blank()
#         )
#     # --- 3. Create Right Plots (Zoomed Clusters) ---
    
#     zoom_plots <- list()
    
#     for (method_name in methods_to_plot) {
        
#         # --- Get Point Data ---
#         if (!method_name %in% names(all_clusterings)) {
#             warning(paste("Method not found in all_clusterings:", method_name))
#             next
#         }
#         pts_df <- all_clusterings[[method_name]]
#         if (is.list(pts_df) && "result_df" %in% names(pts_df)) {
#           pts_df <- pts_df$result_df
#         }
        
#         # --- Get Geometry Data ---
#         if (!method_name %in% names(all_site_geometries)) {
#             warning(paste("Method not found in all_site_geometries:", method_name))
#             next
#         }
#         geom_sf <- all_site_geometries[[method_name]]
#         if (is.null(geom_sf)) {
#             warning(paste("NULL geometry for method:", method_name))
#             next
#         }
        
#         # --- Transform Geometries to WGS84 for plotting ---
#         if (is.na(sf::st_crs(geom_sf))) {
#              sf::st_crs(geom_sf) <- 5070 
#         }
#         geom_sf_wgs84 <- sf::st_transform(geom_sf, crs = wgs84_crs)
#         geom_sf_wgs84$site <- as.factor(geom_sf_wgs84$site)
        
        
#         # +++ FIX START: Robust Intersection Logic +++
        
#         # 1. Create bbox
#         zoom_bbox_sf <- sf::st_bbox(c(
#             xmin = zoom_box$longitude[1], 
#             xmax = zoom_box$longitude[2],
#             ymin = zoom_box$latitude[1], 
#             ymax = zoom_box$latitude[2]
#         ), crs = wgs84_crs)
#         zoom_poly_sfc <- sf::st_as_sfc(zoom_bbox_sf)

#         # 2. Make valid BEFORE intersection (Crucial for fine grids)
#         sf::st_agr(geom_sf_wgs84) = "constant"
#         geom_sf_wgs84 <- sf::st_make_valid(geom_sf_wgs84)

#         # Set agr to constant
#         sf::st_agr(geom_sf_wgs84) <- "constant"

#         # 3. Intersect
#         geom_sf_zoom <- suppressWarnings(
#             sf::st_intersection(geom_sf_wgs84, zoom_poly_sfc)
#         )
        
#         # 4. Post-Intersection Extraction (CRITICAL CHANGE)
#         if (nrow(geom_sf_zoom) > 0) {
            
#             # Ensure we extract POLYGONS from any GEOMETRYCOLLECTIONS created by the cut
#             # This replaces the previous "filter" which was dropping the collections
#             # geom_sf_zoom <- sf::st_collection_extract(geom_sf_zoom, "POLYGON")
            

#             if (inherits(sf::st_geometry(geom_sf_zoom), "sfc_GEOMETRYCOLLECTION") || any(sf::st_geometry_type(geom_sf_zoom) == "GEOMETRYCOLLECTION")) {
#               geom_sf_zoom <- sf::st_collection_extract(geom_sf_zoom, "POLYGON")
#             }
#             # Make valid again after extraction
#             geom_sf_zoom <- sf::st_make_valid(geom_sf_zoom)
            
#             # Filter by area (explicitly stripping units to avoid errors)
#             if (nrow(geom_sf_zoom) > 0) {
#                 # Calculate area
#                 area_vals <- sf::st_area(geom_sf_zoom)
                
#                 # Convert to numeric (strips 'm^2' unit class) for safe comparison
#                 # 1e-6 m^2 is effectively 0, checking for non-empty polygons
#                 geom_sf_zoom <- geom_sf_zoom[as.numeric(area_vals) > 1e-6, ]
#             }
#         }
#         # +++ FIX END +++


#         # --- Filter Point Data ---
#         pts_df$site <- as.factor(pts_df$site)
#         pts_df_zoom <- pts_df[
#             (pts_df$longitude > zoom_box$longitude[1]) &
#             (pts_df$longitude < zoom_box$longitude[2]) &
#             (pts_df$latitude > zoom_box$latitude[1]) &
#             (pts_df$latitude < zoom_box$latitude[2]), ]
        
#         pts_df_distinct <- pts_df_zoom %>%
#             distinct(latitude, longitude, .keep_all = TRUE)
            
#         p_zoom <- ggplot() +
#             geom_raster(
#                 data = base_rast_df, 
#                 aes(x = x, y = y, fill = .data[[elev_col_name]]),
#                 show.legend = FALSE
#             ) +
#             scale_fill_viridis_c(
#                 option = "H",
#                 limits = c(base_rast_min, base_rast_max),
#                 aesthetics = "fill"
#             ) +
#             new_scale_fill() + 
            
#             # --- Site Geometries ---
#             {if (nrow(geom_sf_zoom) > 0) 
#                 geom_sf(
#                     data = geom_sf_zoom, 
#                     aes(fill = site), 
#                     alpha = 0.4,      
#                     color = "black",  
#                     linewidth = 0.25,
#                     show.legend = FALSE,
#                     inherit.aes = FALSE 
#                 )
#             } +
            
#             # --- Clustered points ---
#             geom_point(
#                 data = pts_df_distinct,
#                 aes(x = longitude, y = latitude, fill = site),
#                 shape = 21, size = 2.0, 
#                 color = "black",
#                 show.legend = FALSE
#             ) +
#             scale_fill_discrete() + 
#             coord_sf(
#                 xlim = c(zoom_box$longitude[1], zoom_box$longitude[2]),
#                 ylim = c(zoom_box$latitude[1], zoom_box$latitude[2]),
#                 crs = wgs84_crs,
#                 expand = FALSE 
#             ) +
#             theme_bw() +
#             labs(title = method_name) +
#             theme(
#                 axis.title = element_blank(),
#                 axis.text = element_blank(),
#                 axis.ticks = element_blank(),
#                 plot.title = element_text(size = 10, hjust = 0.5),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank()
#             )
            
#         zoom_plots[[method_name]] <- p_zoom
#     }
    
#     # --- 4. Assemble the Final Plot ---
#     plot_clust <- patchwork::wrap_plots(zoom_plots, ncol = 6)
#     final_plot <- obs_plot + plot_clust +
#         plot_layout(nrow = 1, widths = c(1, 4)) 

#     ggsave(
#         output_path,
#         plot = final_plot,
#         width = 14,
#         height = 8, # 8 for 3 row, 6 for 2 row
#         dpi = 300
#     )
    
#     cat(sprintf("--- Site cluster plot saved to %s ---\n", output_path))
# }

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
    elevation_raster, # INPUT IS NOW ALBERS (res_m)
    methods_to_plot,
    boundary_shp_path,
    output_path,
    zoom_box = list(
        longitude = c(-123.025, -122.992),
        latitude = c(44.085, 44.118)
    )
) {

    # --- 0. Setup Coordinate Systems ---
    albers_crs <- terra::crs(elevation_raster) # Get Albers CRS from input
    wgs84_crs <- sf::st_crs(4326)

    # --- 1. Prepare Rasters ---
    
    # A. Prepare ALBERS Raster (For Right/Zoom Plots - Keeps res_m)
    valid_boundary <- terra::vect(boundary_shp_path)
    valid_boundary_albers <- terra::project(valid_boundary, albers_crs)
    
    region_albers <- terra::crop(elevation_raster[[1]], valid_boundary_albers, mask = TRUE)
    base_rast_df_albers <- as.data.frame(region_albers, xy = TRUE)
    elev_col_name <- names(base_rast_df_albers)[3]
    
    # Get min/max from Albers raster to ensure consistent coloring across both plots
    base_rast_min <- min(base_rast_df_albers[[elev_col_name]], na.rm = TRUE)
    base_rast_max <- max(base_rast_df_albers[[elev_col_name]], na.rm = TRUE)

    # B. Prepare WGS84 Raster (For Left/Overview Plot)
    # Project the Albers region back to WGS84 to align with base_train_df (lat/long)
    region_wgs84 <- terra::project(region_albers, wgs84_crs)
    base_rast_df_wgs84 <- as.data.frame(region_wgs84, xy = TRUE)
    bbox_full <- terra::ext(region_wgs84) # Extent for the left plot

    
    # --- 2. Create Left Plot (Observations in WGS84) ---
    # This matches the "old" style using the WGS84 projected raster
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
        # Annotations (Coordinates are approx WGS84)
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
                hjust = 0.5,
                face = "bold",
                vjust = -5,
                margin = margin(b = -10)
            ),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.margin = margin(t = -250),
            legend.box.margin = margin(0, 0, 0, 0),
            legend.key.width = unit(0.75, "cm"),
            legend.key.height = unit(0.5, "cm"),
            legend.title = element_text(size = 10, vjust = 1),
            legend.text = element_text(size = 8),
            axis.title = element_text(size = 10)
        )

    # --- 3. Create Right Plots (Zoomed Clusters in ALBERS) ---
    
    # Transform Zoom Box to Albers Polygon for cropping
    zoom_poly_wgs84 <- sf::st_as_sfc(sf::st_bbox(c(
        xmin = zoom_box$longitude[1], xmax = zoom_box$longitude[2],
        ymin = zoom_box$latitude[1], ymax = zoom_box$latitude[2]
    ), crs = wgs84_crs))
    
    zoom_poly_albers <- sf::st_transform(zoom_poly_wgs84, albers_crs)
    zoom_bbox_albers <- sf::st_bbox(zoom_poly_albers) # For coord_sf limits

    zoom_plots <- list()

    for (method_name in methods_to_plot) {

        # --- Get Point Data & Transform to Albers ---
        if (!method_name %in% names(all_clusterings)) next
        pts_df <- all_clusterings[[method_name]]
        if (is.list(pts_df) && "result_df" %in% names(pts_df)) pts_df <- pts_df$result_df
        
        # Project points to Albers
        pts_sf <- sf::st_as_sf(pts_df, coords = c("longitude", "latitude"), crs = wgs84_crs)
        pts_sf_albers <- sf::st_transform(pts_sf, albers_crs)
        
        # Extract coordinates for plotting
        pts_coords <- sf::st_coordinates(pts_sf_albers)
        pts_df_albers <- pts_df
        pts_df_albers$x <- pts_coords[,1]
        pts_df_albers$y <- pts_coords[,2]
        pts_df_albers$site <- as.factor(pts_df_albers$site)

        # Filter points to zoom box (using Albers coords)
        pts_df_zoom <- pts_df_albers[
            (pts_df_albers$x > zoom_bbox_albers["xmin"]) &
            (pts_df_albers$x < zoom_bbox_albers["xmax"]) &
            (pts_df_albers$y > zoom_bbox_albers["ymin"]) &
            (pts_df_albers$y < zoom_bbox_albers["ymax"]), 
        ]
        
        # --- Get Geometry Data (Already Albers) ---
        if (!method_name %in% names(all_site_geometries)) next
        geom_sf <- all_site_geometries[[method_name]]
        if (is.null(geom_sf)) next
        
        # Ensure geom_sf has CRS set
        if (is.na(sf::st_crs(geom_sf))) sf::st_crs(geom_sf) <- albers_crs
        geom_sf$site <- as.factor(geom_sf$site)

        # Crop geometries to Albers Zoom Box
        # Use st_crop or intersection. st_crop is faster for bboxes.
        geom_sf_zoom <- suppressWarnings(sf::st_crop(geom_sf, zoom_bbox_albers))
        
        # Clean up collections if generated
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
            scale_fill_discrete() +
            
            # Enforce Albers Coordinates
            coord_sf(
                xlim = c(zoom_bbox_albers["xmin"], zoom_bbox_albers["xmax"]),
                ylim = c(zoom_bbox_albers["ymin"], zoom_bbox_albers["ymax"]),
                crs = albers_crs,
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
        height = 8,
        dpi = 300
    )

    cat(sprintf("--- Site cluster plot saved to %s ---\n", output_path))
}
