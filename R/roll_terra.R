#' In-Memory Rolling Time Series Averaging (Internal)
#' @keywords internal
roll_terra <- function(x, window_size, freq, step, fun, output_names_list, output_dir, gdal_opt, overwrite) {
  
  total_cycles <- nlyr(x) / freq
  start_units <- seq(1, total_cycles - window_size + 1, by = step)
  
  all_results_stack <- list()
  
  # Loop over each window
  for (i in seq_along(start_units)) {
    start_y <- start_units[i]
    end_y <- start_y + window_size - 1
    
    rlang::inform(glue::glue("Processing window: Cycle {start_y} to {end_y}"))
    
    # Get layers for the current window
    start_layer_idx <- ((start_y - 1) * freq) + 1
    end_layer_idx <- end_y * freq
    sub_stack <- x[[start_layer_idx:end_layer_idx]]
    
    # Create index for tapp (e.g., 1,2,...,12, 1,2,...,12, ...)
    grouping_idx <- rep(1:freq, times = window_size)
    
    # Calculate the average for this window
    avg_rast_window <- terra::tapp(sub_stack, index = grouping_idx, fun = fun, na.rm = TRUE)
    
    # Use the pre-generated names for this window's layers
    names(avg_rast_window) <- output_names_list[[i]]
    
    all_results_stack[[i]] <- avg_rast_window
  }
  
  # Combine all results into a single SpatRaster
  final_stack <- terra::rast(all_results_stack)
  
  output_files <- file.path(output_dir, paste0(names(final_stack), ".tif"))
  
  rlang::inform("Writing final GeoTIFFs...")
  terra::writeRaster(
    final_stack,
    filename = output_files,
    overwrite = overwrite,
    gdal = gdal_opt
  )
  
  return(output_files)
}