#' In-Memory Time Series Averaging (Internal)
#' @keywords internal
average_terra <- function(x, index, output_names, output_dir, gdal_opt, overwrite) {
  rlang::inform("Calculating averages using terra::tapp...")
  
  avg_rast <- terra::tapp(x, index = index, fun = "mean", na.rm = TRUE)
  names(avg_rast) <- output_names
  
  output_files <- file.path(output_dir, paste0(names(avg_rast), ".tif"))
  
  rlang::inform("Writing final GeoTIFFs...")
  terra::writeRaster(
    avg_rast,
    filename = output_files,
    overwrite = overwrite,
    gdal = gdal_opt
  )
  return(output_files)
}