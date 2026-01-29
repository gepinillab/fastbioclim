#' In-Memory Time Series Aggregation (Internal)
#' @keywords internal
aggregate_terra <- function(x, index, output_names, output_dir, gdal_opt, overwrite, verbose, aggregation_type = "mean") {
  
  # --- Input Validation ---
  if (!aggregation_type %in% c("mean", "sum")) {
    stop("Invalid 'aggregation_type'. Must be one of 'mean' or 'sum'.")
  }
  if (verbose) message(glue::glue("Calculating aggregations ('{aggregation_type}') using terra::tapp..."))
  agg_rast <- terra::tapp(x, index = index, fun = aggregation_type, na.rm = TRUE)
  names(agg_rast) <- output_names
  output_files <- file.path(output_dir, paste0(names(agg_rast), ".tif"))
  if (verbose) message("Writing final GeoTIFFs...")
  terra::writeRaster(
    agg_rast,
    filename = output_files,
    overwrite = overwrite,
    gdal = gdal_opt
  )
  return(output_files)
}