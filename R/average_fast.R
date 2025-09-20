#' Tiled, Out-of-Core Time Series Averaging (Internal)
#' @keywords internal
average_fast <- function(paths, index, output_names, user_region, tile_degrees, output_dir) {
  
  # --- 1. Setup Environment ---
  qs_dir <- file.path(output_dir, paste0("avg_qs_", basename(tempfile(pattern = ""))))
  dir.create(qs_dir, recursive = TRUE)
  # --- 2. Define Processing Region and Replicate Full Geometry Logic ---
  # --- Check Geometry Consistency ---
  message("Checking geometry of input rasters...")
  tryCatch({
    ref_rast <- terra::rast(paths[1])
    ref_crs <- terra::crs(ref_rast)
    ref_ext <- terra::ext(ref_rast)
    if (length(paths) > 1) {
      for(i in 2:length(paths)) {
        current_rast_info <- try(terra::rast(paths[i]), silent = TRUE)
        if (inherits(current_rast_info, "try-error")) {
          stop("Failed to read header for: ", paths[i])
        }
        terra::compareGeom(ref_rast, current_rast_info, stopOnError = TRUE, messages = FALSE)
      }
    }
  }, error = function(e) {
    stop("Input rasters (including static indices) do not have the same geometry. ", e$message)
  })
  rm(ref_rast)

  # --- 3. Define Processing Region ---
  base_map <- NULL
  if (!is.null(user_region)) {
    message("Using user-provided region.")
    if (inherits(user_region, "SpatVector")) base_map <- sf::st_as_sf(user_region)
    else if (inherits(user_region, "sf") || inherits(user_region, "sfc")) base_map <- sf::st_as_sf(sf::st_geometry(user_region))
    else stop("'user_region' must be an sf object or a terra SpatVector.")
    if (!all(sf::st_is_valid(base_map))) base_map <- sf::st_make_valid(base_map)
    # Transform user_region CRS to match raster CRS if needed
    if(sf::st_crs(base_map) != ref_crs) {
      base_map <- sf::st_transform(base_map, ref_crs)
    }
    # Check if user region overlaps with raster extent
    region_bbox <- sf::st_bbox(base_map)
    raster_bbox <- c(xmin = ref_ext[1], ymin = ref_ext[3], xmax = ref_ext[2], ymax = ref_ext[4])
    # Simple bounding box overlap check
    if (region_bbox$xmax < raster_bbox["xmin.xmin"] || region_bbox$xmin > raster_bbox["xmax.xmax"] ||
        region_bbox$ymax < raster_bbox["ymin.ymin"] || region_bbox$ymin > raster_bbox["ymax.ymax"]) {
          stop("Provided user_region does not overlap with the extent of the input rasters.")
      }
        
    } else {
      message("No user_region provided. Using the full extent of input rasters.")
      # Create an sf polygon from the raster extent obtained earlier
      base_map <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(ref_ext), crs = ref_crs))
    }
        
    # --- Extract and Save Template Geometry Info ---
    ref_rast_geom <- terra::rast(paths[1])
    original_extent_vec <- as.vector(terra::ext(ref_rast_geom))
    original_dims_vec <- c(terra::nrow(ref_rast_geom), terra::ncol(ref_rast_geom))
    original_crs_txt <- terra::crs(ref_rast_geom, proj = TRUE)
    original_ncol <- original_dims_vec[2]
    original_res <- terra::res(ref_rast_geom)

    # Determine target geometry
    target_extent_vec <- original_extent_vec
    target_dims_vec <- original_dims_vec
    target_crs_txt <- original_crs_txt
    target_ncol <- original_ncol
    target_res <- original_res

    if (!is.null(user_region)) {
      target_template_rast <- NULL
      tryCatch({
        target_template_rast <- terra::crop(ref_rast_geom, terra::vect(base_map), mask = TRUE)
        if (terra::ncell(target_template_rast) > 0) {
          target_extent_vec <- as.vector(terra::ext(target_template_rast))
          target_dims_vec <- c(terra::nrow(target_template_rast), terra::ncol(target_template_rast))
          target_ncol <- target_dims_vec[2]
          target_res <- terra::res(target_template_rast)
        } else {
          warning("Cropping resulted in empty raster, using original geometry.")
        }
        rm(target_template_rast)
      }, error = function(e){
        warning("Could not derive exact cropped geometry, using full raster extent info. Error: ", e$message)
      })
    }

    template_info <- list(
      original_geom = list(
        extent = original_extent_vec,
        dimensions = original_dims_vec,
        crs = original_crs_txt
      ),
      target_geom = list(
        extent = target_extent_vec,
        dimensions = target_dims_vec,
        crs = target_crs_txt,
        res = target_res
      )
    )
    # --- Calculate Translation Parameters ---
    if (!is.null(user_region)) {
    # Use the original raster *header* info (ref_rast_geom) for offset calculation
    target_xmin <- template_info$target_geom$extent[1] + template_info$target_geom$res[1] / 2
    target_ymax <- template_info$target_geom$extent[4] - template_info$target_geom$res[2] / 2
    # Handle potential floating point inaccuracies by finding the *nearest* cell
    col_offset <- tryCatch(
      terra::colFromX(ref_rast_geom, target_xmin) - 1L,
      error = function(e) {
        warning("Could not precisely determine column offset using terra::colFromX, possibly due to edge alignment. Assuming 0 offset. Error: ", e$message)
        0L # Default to 0 if calculation fails
      }
    )
    row_offset <- tryCatch(
      terra::rowFromY(ref_rast_geom, target_ymax) - 1L,
      error = function(e) {
        warning("Could not precisely determine row offset using terra::rowFromY, possibly due to edge alignment. Assuming 0 offset. Error: ", e$message)
        0L # Default to 0 if calculation fails
      }
    )
  
    # Ensure offsets are non-negative (should be if target is within original)
    col_offset <- max(0L, col_offset)
    row_offset <- max(0L, row_offset)
  
    # --- Create the specific translation function ---
    translate_cell <- define_translate(
      ncol_src = original_ncol,
      ncol_tgt = target_ncol,
      row_offset = row_offset,
      col_offset = col_offset
    )
  }
    
  # Save the template info within that directory
  template_info_file <- file.path(qs_dir, "template_info.qs2")
  tryCatch({
    qs2::qs_save(template_info, template_info_file)
  }, error = function(e){
    stop("Failed to save template geometry information: ", e$message)
  })
  rm(ref_rast_geom)
  
  
  # --- 3. Create Spatial Tiles ---
  sf::sf_use_s2(FALSE)
  grid_bbox <- sf::st_bbox(base_map)
  rtt_grid <- sf::st_make_grid(grid_bbox, cellsize = tile_degrees, 
                               what = "polygons", crs = sf::st_crs(base_map))
  rtt <- sf::st_intersection(base_map, rtt_grid)
  rtt <- sf::st_collection_extract(rtt, c("POLYGON"))
  rtt <- sf::st_as_sf(rtt)
  rtt <- rtt[!sf::st_is_empty(rtt),]
  ntiles <- nrow(rtt)
  if (ntiles == 0) stop("No overlapping tiles found for the processing area.")
  message("Rasters divided into ", ntiles, " tiles for processing.")
  
  # --- 4. Parallel Processing ---
  p <- progressr::progressor(steps = ntiles)
  unique_indices <- sort(unique(index))
  n_units_out <- length(unique_indices)
  export_vars <- c("paths", "path_variables", "rtt", "ntiles", "unique_indices", 
                   "n_units_out", "output_names", "qs_dir", "translate_cell")
  future.apply::future_lapply(seq_len(ntiles), function(i) {
    p(message = sprintf("Processing tile %d/%d", i, ntiles))
    tile_geom <- sf::st_geometry(rtt[i, ])
    
    evars_stack_tile <- tryCatch({ terra::rast(paths) }, error = function(e) { NULL })
        if (is.null(evars_stack_tile)) { 
          warning(sprintf("Tile %d: Failed load", i))
          return(NULL)
        }

    extracted <- tryCatch(
      exactextractr::exact_extract(evars_stack_tile, tile_geom, include_cell = TRUE),
      error = function(e) NULL
    )
    if (is.null(extracted) || nrow(extracted[[1]]) == 0) return(NULL)
    df <- extracted[[1]][!is.na(extracted[[1]][[names(extracted[[1]])[1]]]), ]
    if(nrow(df) == 0) return(NULL)
    nonaID <- which(!is.na(df[[2]]))
    if (length(nonaID) == 0L) return(NULL)
    cell_ids_source <- df$cell[nonaID]
    if(!is.null(user_region)) {
      cell_ids_source <- translate_cell(cell_ids_source)
      noMap <- !is.na(cell_ids_source)
      cell_ids_source <- cell_ids_source[noMap]
    }
    value_cols <- names(df)[!names(df) %in% c("cell", "coverage_fraction")]
    pixel_matrix <- Rfast::data.frame.to_matrix(df[nonaID[noMap], value_cols, drop = FALSE])
    # Calculate average for each output group and save it to a separate .qs2 file
    for (j in 1:n_units_out) {
      current_idx_val <- unique_indices[j]
      current_output_name <- output_names[j]
      
      cols_to_avg <- which(index == current_idx_val)
      avg_values <- Rfast::rowmeans(pixel_matrix[, cols_to_avg, drop = FALSE])
      
      # Create 2-column data frame that write_layers expects: value, cell
      tile_result <- data.frame(value = avg_values, cell = cell_ids_source)
      # Filename format: [variable_name]_[tile_number].qs2
      qs_filename <- file.path(qs_dir, paste0(current_output_name, "_", i, ".qs2"))
      qs2::qs_save(tile_result, qs_filename)
    }
    return(NULL)
  }, 
  future.seed = TRUE, future.globals = export_vars,
  future.packages = c("terra", "exactextractr", "Rfast", "qs", "sf"))
  
  rlang::inform("Tiled computation finished. Assembly will be handled by write_layers.")
  return(qs_dir)
}