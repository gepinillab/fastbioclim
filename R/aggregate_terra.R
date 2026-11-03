#' In-Memory Time Series Aggregation (Internal)
#' @param x A `terra::SpatRaster` with the time series to aggregate.
#' @param index Vector giving the group each layer belongs to
#'   (see `terra::tapp()`).
#' @param output_names Character vector with the names of the aggregated
#'   output layers.
#' @param output_dir Character, path to save final rasters.
#' @param gdal_opt Character vector of GDAL options for writing.
#' @param overwrite Logical, whether to overwrite existing files.
#' @param verbose Logical, If `TRUE`, prints messages.
#' @param aggregation_type Character, either `"mean"` or `"sum"`.
#' @keywords internal
aggregate_terra <- function(x, index, output_names, output_dir, gdal_opt, overwrite, verbose, aggregation_type = "mean") {
  
  # --- Input Validation ---
  if (!aggregation_type %in% c("mean", "sum")) {
    stop("Invalid 'aggregation_type'. Must be one of 'mean' or 'sum'.")
  }
  if (verbose) message(glue::glue("Calculating aggregations ('{aggregation_type}') using terra::tapp..."))
  # tapp() orders groups by first appearance unless index is a factor; use the
  # sorted order that output_names (and the tiled backend) follow.
  index <- factor(index, levels = sort(unique(index)))
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