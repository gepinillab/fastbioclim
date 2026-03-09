#' Calculate Sums for SpatRasters
#'
#' Calculates temporal sums for a multi-layer SpatRaster.
#' This function serves as a smart wrapper, automatically selecting between
#' an in-memory (`terra`) or out-of-core (`tiled`) workflow based on data size.
#'
#' @param x A `terra::SpatRaster` object with multiple layers representing a time series.
#' @param index A numeric or integer vector defining the grouping for aggregation.
#'   Its length must equal the number of layers in `x`. For example, to sum
#'   365 daily layers into 12 monthly totals, `index` would group the days by month.
#' @param output_names A character vector of names for the output layers. Its
#'   length must equal the number of unique groups in `index`. If `NULL`, names
#'   like "sum_unit_1" are generated.
#' @param output_dir The directory where the final summed raster layers
#'   will be saved as GeoTIFF files.
#' @param user_region (Optional) An `sf` or `terra::SpatVector` object. If provided,
#'   the input raster `x` is clipped and masked to this region before processing.
#'   The output raster's extent is the same of the `user_region`.
#' @param method The processing method: "auto", "tiled", or "terra".
#' @param tile_degrees (Tiled method only) The approximate size of processing tiles.
#' @param gdal_opt (Optional) GDAL creation options for the output GeoTIFFs.
#' @param overwrite Logical. If `FALSE` (default), stops if output files exist.
#' @param verbose Logical, If `TRUE`, prints messages.
#' @return A `terra::SpatRaster` object pointing to the newly created files, with the following characteristics:
#'   \itemize{
#'     \item **Number of layers:** The number of layers will be equal to the number of unique values in the `index` argument.
#'     \item **Layer names:** Layer names are determined by the `output_names` argument. If `NULL`, they will be generated automatically (e.g., 'sum_unit_01', 'sum_unit_02', etc.).
#'     \item **Extent:** If `user_region` is provided, the extent of the output raster will be clipped to match that region. Otherwise, the extent will be the same as the input raster `x`.
#'   }
#' @examples
#' \donttest{
#' # The example raster "prcp.tif" is included in the package's `inst/extdata` directory.
#' # Load example data from Lesotho (Monthly time series from 2016-01 to 2020-12)
#' raster_path <- system.file("extdata", "prcp.tif", package = "fastbioclim")
#' # Load the SpatRaster from the file
#' prcp_ts <- terra::rast(raster_path)
#' # To calculate the total annual precipitation, we group the 60 monthly layers by year.
#' annual_index <- rep(2016:2020, each = 12)
#' # Define a temporary directory for the output files
#' output_dir <- file.path(tempdir(), "annual_prcp_sum")
#' # Run the calculate_sum function
#' annual_sum <- calculate_sum(
#'   x = prcp_ts,
#'   index = annual_index,
#'   output_names = "prcp_sum",
#'   output_dir = output_dir,
#'   overwrite = TRUE,
#'   verbose = FALSE
#' )
#' # Print the resulting SpatRaster summary (should have 5 layers)
#' print(annual_sum)
#' # Clean up the created files
#' unlink(output_dir, recursive = TRUE)
#' }
#' @export
calculate_sum <- function(x,
                          index,
                          output_names = NULL,
                          output_dir = tempdir(),
                          user_region = NULL,
                          method = c("auto", "tiled", "terra"),
                          tile_degrees = 5,
                          gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
                          overwrite = FALSE,
                          verbose = TRUE) {

  # --- 1. Argument Validation and Setup ---
  method <- match.arg(method)
  if (!inherits(x, "SpatRaster")) stop("Input 'x' must be a SpatRaster.")
  n_in <- nlyr(x)

  # Validate index
  if (missing(index) || !is.numeric(index)) stop("A numeric 'index' vector is required.")
  if (length(index) != n_in) {
    stop(glue::glue("Length of 'index' ({length(index)}) must match number of layers in 'x' ({n_in})."))
  }

  unique_groups <- sort(unique(index))
  n_units_out <- length(unique_groups)

  # Validate or create output_names
  if (is.null(output_names)) {
    suffix_name <- "sum_unit_"
    output_names <- paste0("sum_unit_", sprintf("%02d", unique_groups))
  } else {
    suffix_name <- output_names
    output_names <- paste0(output_names, "_", sprintf("%02d", unique_groups))
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --- 2. Fail-Fast Pre-Check for Existing Files ---
  if (!overwrite) {
    expected_filepaths <- file.path(output_dir, paste0(output_names, ".tif"))
    if (any(file.exists(expected_filepaths))) {
      stop("Output files already exist and `overwrite` is FALSE.")
    }
  }

  # --- 3. The "Auto" Decision Logic ---
  use_terra_workflow <- FALSE
  if (method == "terra") {
    use_terra_workflow <- TRUE
    if (verbose) message("User forced 'terra' (in-memory) workflow.")
  } else if (method == "tiled") {
    use_terra_workflow <- FALSE
    if (verbose) message("User forced 'tiled' workflow.")
  } else { # "auto"
    if (verbose) message("Using 'auto' method to select workflow...")
    template_rast <- if (is.null(user_region)) x[[1]] else crop(x[[1]], user_region)
    mem_info <- terra::mem_info(template_rast, n = nlyr(x) + n_units_out, print = FALSE)
    if (mem_info["fits_mem"] == 1) {
      use_terra_workflow <- TRUE
      if (verbose) message("Data appears to fit in memory. Selecting 'terra' workflow.")
    } else {
      use_terra_workflow <- FALSE
      if (verbose) message("Data is too large for memory. Selecting 'tiled' workflow.")
    }
  }

  # --- 4. Delegate to the Correct Workflow ---
  if (use_terra_workflow) {
    final_raster <- if (!is.null(user_region)) {
      if (verbose) message("Performing crop operation for 'terra' workflow...")
      terra::crop(x, user_region, mask = TRUE)
    } else {
      x
    }

    result_paths <- aggregate_terra(
      x = final_raster,
      index = index,
      output_names = output_names,
      output_dir = output_dir,
      gdal_opt = gdal_opt,
      overwrite = overwrite,
      verbose = verbose,
      aggregation_type = "sum"
    )

  } else { # Tiled workflow
    if (all(terra::inMemory(x))) {
      stop("The 'tiled' workflow requires the input SpatRaster to point to a file on disk.")
    }

    temp_qs_dir <- aggregate_fast(
      paths = terra::sources(x),
      index = index,
      output_names = output_names,
      user_region = user_region,
      tile_degrees = tile_degrees,
      output_dir = output_dir,
      verbose = verbose,
      aggregation_type = "sum"
    )

    # Use write_layers to assemble the results
    result_paths <- write_layers(
      input_dir = temp_qs_dir,
      output_dir = output_dir,
      file_pattern = suffix_name,
      gdal_opt = gdal_opt,
      overwrite = overwrite,
      verbose = verbose
    )
  }
  if (verbose) message(paste("Processing complete. Final rasters are in:", normalizePath(output_dir)))
  return(terra::rast(result_paths))
}