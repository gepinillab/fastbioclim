#' Tiled, Out-of-Core Rolling Time Series Averaging (Internal) - CORRECTED
#' @keywords internal
roll_fast <- function(paths, window_size, freq, step, fun, output_names_list, user_region, tile_degrees, output_dir) {
  
  # --- 1. Setup Environment ---
  qs_dir <- file.path(output_dir, paste0("roll_avg_qs_", basename(tempfile(pattern = ""))))
  dir.create(qs_dir, recursive = TRUE)

  # --- 2. Geometry, Region, and Tiling Setup ---
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
  base_map <- NULL
  if (!is.null(user_region)) {
    message("Using user-provided region.")
    if (inherits(user_region, "SpatVector")) base_map <- sf::st_as_sf(user_region)
    else if (inherits(user_region, "sf") || inherits(user_region, "sfc")) base_map <- sf::st_as_sf(sf::st_geometry(user_region))
    else stop("'user_region' must be an sf object or a terra SpatVector.")
    if (!all(sf::st_is_valid(base_map))) base_map <- sf::st_make_valid(base_map)
    if(sf::st_crs(base_map) != ref_crs) {
      base_map <- sf::st_transform(base_map, ref_crs)
    }
    region_bbox <- sf::st_bbox(base_map)
    raster_bbox <- c(xmin = ref_ext[1], ymin = ref_ext[3], xmax = ref_ext[2], ymax = ref_ext[4])
    if (region_bbox$xmax < raster_bbox["xmin.xmin"] || region_bbox$xmin > raster_bbox["xmax.xmax"] ||
        region_bbox$ymax < raster_bbox["ymin.ymin"] || region_bbox$ymin > raster_bbox["ymax.ymax"]) {
          stop("Provided user_region does not overlap with the extent of the input rasters.")
      }
    } else {
      message("No user_region provided. Using the full extent of input rasters.")
      base_map <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(ref_ext), crs = ref_crs))
    }
    ref_rast_geom <- terra::rast(paths[1])
    original_extent_vec <- as.vector(terra::ext(ref_rast_geom))
    original_dims_vec <- c(terra::nrow(ref_rast_geom), terra::ncol(ref_rast_geom))
    original_crs_txt <- terra::crs(ref_rast_geom, proj = TRUE)
    original_ncol <- original_dims_vec[2]
    original_res <- terra::res(ref_rast_geom)
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
    if (!is.null(user_region)) {
    target_xmin <- template_info$target_geom$extent[1] + template_info$target_geom$res[1] / 2
    target_ymax <- template_info$target_geom$extent[4] - template_info$target_geom$res[2] / 2
    col_offset <- tryCatch(
      terra::colFromX(ref_rast_geom, target_xmin) - 1L,
      error = function(e) {
        warning("Could not precisely determine column offset using terra::colFromX, possibly due to edge alignment. Assuming 0 offset. Error: ", e$message)
        0L
      }
    )
    row_offset <- tryCatch(
      terra::rowFromY(ref_rast_geom, target_ymax) - 1L,
      error = function(e) {
        warning("Could not precisely determine row offset using terra::rowFromY, possibly due to edge alignment. Assuming 0 offset. Error: ", e$message)
        0L
      }
    )
    col_offset <- max(0L, col_offset)
    row_offset <- max(0L, row_offset)
    translate_cell <- define_translate(
      ncol_src = original_ncol,
      ncol_tgt = target_ncol,
      row_offset = row_offset,
      col_offset = col_offset
    )
  }
  template_info_file <- file.path(qs_dir, "template_info.qs")
  tryCatch({
    qs::qsave(template_info, template_info_file)
  }, error = function(e){
    stop("Failed to save template geometry information: ", e$message)
  })
  rm(ref_rast_geom)
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
  
  total_cycles <- length(paths) / freq
  start_units <- seq(1, total_cycles - window_size + 1, by = step)

  export_vars <- c("paths", "rtt", "ntiles", "window_size", "freq", "fun",
                   "start_units", "output_names_list", "qs_dir", "translate_cell", "user_region")

  future.apply::future_lapply(seq_len(ntiles), function(i) {
    p(message = sprintf("Processing tile %d/%d", i, ntiles))
    tile_geom <- sf::st_geometry(rtt[i, ])

    # Lógica de extracción idéntica a la tuya
    evars_stack_tile <- tryCatch({ terra::rast(paths) }, error = function(e) NULL)
    if (is.null(evars_stack_tile)) return(NULL)
    extracted <- tryCatch(
      exactextractr::exact_extract(evars_stack_tile, tile_geom, include_cell = TRUE),
      error = function(e) NULL
    )
    if (is.null(extracted) || nrow(extracted[[1]]) == 0) return(NULL)
    df <- extracted[[1]][!is.na(extracted[[1]][[names(extracted[[1]])[1]]]), ]
    if(nrow(df) == 0) return(NULL)
    nonaID <- which(!is.na(df[[2]]))
    if (length(nonaID) == 0L) return(NULL)
    cell_ids_source_orig <- df$cell[nonaID]
    if(!is.null(user_region)) {
      cell_ids_target <- translate_cell(cell_ids_source_orig)
      noMap <- !is.na(cell_ids_target)
      cell_ids_source <- cell_ids_target[noMap]
      value_cols <- names(df)[!names(df) %in% c("cell", "coverage_fraction")]
      pixel_matrix <- Rfast::data.frame.to_matrix(df[nonaID[noMap], value_cols, drop = FALSE])
    } else {
      cell_ids_source <- cell_ids_source_orig
      value_cols <- names(df)[!names(df) %in% c("cell", "coverage_fraction")]
      pixel_matrix <- Rfast::data.frame.to_matrix(df[nonaID, value_cols, drop = FALSE])
    }
    
    # *** LÓGICA CORREGIDA PARA EL PROMEDIO MÓVIL ***
    for (win_idx in seq_along(start_units)) {
      start_y <- start_units[win_idx]
      end_y <- start_y + window_size - 1
      
      start_col_idx <- ((start_y - 1) * freq) + 1
      end_col_idx <- end_y * freq
      
      window_pixel_matrix <- pixel_matrix[, start_col_idx:end_col_idx, drop = FALSE]
      
      # *** INICIO DEL BLOQUE CORREGIDO ***
      # Usaremos un bucle (pequeño, de 1 a freq) y Rfast::rowmeans/rowsums
      
      period_results <- vector("list", freq)
      
      for (p_idx in 1:freq) {
        # Obtener las columnas para este período (ej. todos los eneros de la ventana)
        period_cols <- seq(from = p_idx, to = ncol(window_pixel_matrix), by = freq)
        
        # Subconjunto de la matriz de la ventana para este período
        period_matrix <- window_pixel_matrix[, period_cols, drop = FALSE]
        
        # Calcular la estadística usando una función Rfast existente y rápida
        if (fun == "mean") {
          period_results[[p_idx]] <- Rfast::rowmeans(period_matrix)
        } else if (fun == "sum") {
          period_results[[p_idx]] <- Rfast::rowsums(period_matrix)
        } else {
          # Fallback para otras funciones (un poco más lento)
          period_results[[p_idx]] <- apply(period_matrix, 1, FUN = fun, na.rm = TRUE)
        }
      }
      
      # Combinar los vectores de resultados en una matriz
      avg_matrix_for_window <- do.call(cbind, period_results)
      
      # *** FIN DEL BLOQUE CORREGIDO ***

      output_names_for_window <- output_names_list[[win_idx]]
      
      for (p_idx in 1:freq) {
        current_output_name <- output_names_for_window[p_idx]
        avg_values <- avg_matrix_for_window[, p_idx]
        
        tile_result <- data.frame(value = avg_values, cell = cell_ids_source)
        qs_filename <- file.path(qs_dir, paste0(current_output_name, "_", i, ".qs"))
        qs::qsave(tile_result, qs_filename)
      }
    }
    return(NULL)
  }, 
  future.seed = TRUE, future.globals = export_vars,
  future.packages = c("terra", "exactextractr", "Rfast", "qs", "sf", "rlang", "glue"))
  
  rlang::inform("Tiled computation finished. Assembly will be handled by write_layers.")
  return(qs_dir)
}