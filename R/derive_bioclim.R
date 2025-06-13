#' Derive Comprehensive Bioclimatic Variables
#'
#' Calculates up to 35 bioclimatic variables from average monthly climate SpatRasters
#' (or other temporal units). This function serves as a smart wrapper that
#' automatically selects the most efficient processing workflow (in-memory vs.
#' tiled) based on data size and user-defined region of interest.
#'
#' @details
#' This function unifies two processing backends. The `method` argument controls
#' which is used:
#' \itemize{
#'   \item `"auto"`: (Default) Intelligently chooses between "terra" and "tiled"
#'     based on estimated memory requirements.
#'   \item `"terra"`: Forces the fast, in-memory workflow. May fail on very large datasets.
#'   \item `"tiled"`: Forces the memory-safe, out-of-core workflow. Ideal for
#'     very large datasets. Requires that the input SpatRasters point to files on disk.
#' }
#' Period-based variables (e.g., BIO8, BIO10) are calculated using a moving
#' window defined by `period_length`.
#'
#' @section Static Indices:
#' For advanced use cases, such as time-series analysis or defining specific
#' seasons, you can provide pre-calculated index rasters to override the
#' dynamic calculations. These are passed as named `SpatRaster` objects via the
#' `...` argument (e.g., `warmest_period = my_warmest_idx_rast`). The wrapper
#' function automatically handles passing them to the appropriate workflow.
#'
#' When using the "tiled" workflow, these static index rasters **must** be
#' file-backed (i.e., not held entirely in memory). Supported static indices
#' include:
#' \itemize{
#'   \item `warmest_unit`, `coldest_unit`, `wettest_unit`, `driest_unit`
#'   \item `high_rad_unit`, `low_rad_unit`, `high_mois_unit`, `low_mois_unit`
#'   \item `warmest_period`, `coldest_period`, `wettest_period`, `driest_period`
#'   \item `high_mois_period`, `low_mois_period`
#' }
#'
#' @param bios Numeric vector specifying which bioclimatic variables (1-35) to compute.
#' @param tmin,tmax,prcp,tavg,srad,mois (Optional) `terra::SpatRaster` objects
#'   containing the climate data for each temporal unit (e.g., 12 monthly layers).
#'   All provided rasters must have the same geometry and number of layers.
#' @param output_dir The directory where the final bioclimatic variable rasters
#'   will be saved. The directory will be created if it does not exist. The default
#'   is temporal directory created by `tempdir`. 
#' @param period_length Integer. The number of temporal units (e.g., months) that
#'   define a "period" for calculating summary variables like BIO8 (Mean Temp of
#'   Wettest Quarter). Defaults to 3, representing quarters for monthly data.
#' @param circular Logical. If `TRUE` (the default), period calculations will
#'   wrap around the beginning and end of the time series (e.g., for monthly
#'   data, Dec-Jan-Feb is considered a valid period).
#' @param user_region (Optional) An `sf` or `terra::SpatVector` object defining the
#'   area of interest. If provided, all calculations will be clipped to this region.
#' @param method The processing method. See Details for more information.
#' @param tile_degrees (Tiled method only) The approximate size of processing tiles
#'   in degrees. Ignored if the 'terra' workflow is used.
#' @param gdal_opt (Optional) A character vector of GDAL creation options for the
#'   output GeoTIFF files. Controls compression, threading, etc.
#' @param overwrite (Optional) Logical. If `FALSE` (the default), the function will
#'   stop immediately if any target output files already exist.
#' @param ... Additional arguments, primarily for passing static index rasters. See
#'   the "Static Indices" section for details.
#'
#' @return An SpatRaster with 35 bioclimatic variables or a subset of them:
#' \describe{
#'   \item{bio01}{Mean Temperature of Units}
#'   \item{bio02}{Mean Diurnal Range}
#'   \item{bio03}{Isothermality}
#'   \item{bio04}{Temperature Seasonality}
#'   \item{bio05}{Max Temperature of Warmest Unit}
#'   \item{bio06}{Min Temperature of Coldest Unit}
#'   \item{bio07}{Temperature Range of Units}
#'   \item{bio08}{Mean Temperature of Wettest Period}
#'   \item{bio09}{Mean Temperature of Driest Period}
#'   \item{bio10}{Mean Temperature of Warmest Period}
#'   \item{bio11}{Mean Temperature of Coldest Period}
#'   \item{bio12}{Precipitation Sum}
#'   \item{bio13}{Precipitation of Wettest Unit}
#'   \item{bio14}{Precipitation of Driest Unit}
#'   \item{bio15}{Precipitation Seasonality}
#'   \item{bio16}{Precipitation of Wettest Period}
#'   \item{bio17}{Precipitation of Driest Period}
#'   \item{bio18}{Precipitation of Warmest Period}
#'   \item{bio19}{Precipitation of Coldest Period}
#'   \item{bio20}{Mean Radiation of Units}
#'   \item{bio21}{Highest Radiation Unit}
#'   \item{bio22}{Lowest Radiation Unit}
#'   \item{bio23}{Radiation Seasonality}
#'   \item{bio24}{Radiation of Wettest Period}
#'   \item{bio25}{Radiation of Driest Period}
#'   \item{bio26}{Radiation of Warmest Period}
#'   \item{bio27}{Radiation of Coldest Period}
#'   \item{bio28*}{Mean Moisture Content Of Units}
#'   \item{bio29*}{Highest Moisture Content Unit}
#'   \item{bio30*}{Lowest Moisture Content Unit}
#'   \item{bio31*}{Moisture Content Seasonality}
#'   \item{bio32*}{Mean Moisture Content of Most Moist Period}
#'   \item{bio33*}{Mean Moisture Content of Least Moist Period}
#'   \item{bio34*}{Mean Moisture Content of Warmest Period}
#'   \item{bio35*}{Mean Moisture Content of Coldest Period}
#' }
#' @note 
#' *The original moisture variables proposed in the ANUCLIM manual are based 
#' on the Moisture Index (MI). However, this function allows users to calculate 
#' moisture-based bioclimatic variables using other units of moisture 
#' as inputs, offering greater flexibility in input data usage.
#' 
#' @export
#' @seealso `validate_climate_inputs()` to check data integrity before processing.
#' @references
#' O’Donnell, M. S., & Ignizio, D. A. (2012). Bioclimatic predictors for supporting ecological applications in the conterminous United States.
#' ANUCLIM 6.1 User Guide. Centre for Resource and Environmental Studies, The Australian National University.
#' @examples
#' # This is a conceptual example, requires data setup
#' if (FALSE) {
#'   # Assume tmin_rast, tmax_rast, prcp_rast are 12-layer SpatRasters
#'   bioclim_vars <- derive_bioclim(
#'     bios = 1:19,
#'     tmin = tmin_rast,
#'     tmax = tmax_rast,
#'     prcp = prcp_rast,
#'     output_dir = "./bioclim_output",
#'     overwrite = TRUE
#'   )
#'   plot(bioclim_vars[[c("bio01", "bio12")]])
#' }
#'
derive_bioclim <- function(bios,
  tmin = NULL, tmax = NULL, prcp = NULL, tavg = NULL, srad = NULL, mois = NULL,
  output_dir = tempdir(),
  period_length = 3,
  circular = TRUE,
  user_region = NULL,
  method = c("auto", "tiled", "terra"),
  tile_degrees = 5,
  temp_dir = tempdir(),
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE,
  ...) {

  # --- 1. Argument Capture and Initial Setup ---
  method <- match.arg(method)
  dot_args <- list(...)

  # --- 2. Fail-Fast Pre-Check for Existing Files ---
  if (!overwrite) {
    expected_filenames <- paste0("bio", sprintf("%02d", sort(unique(bios))), ".tif")
    expected_filepaths <- file.path(output_dir, expected_filenames)
    existing_files <- expected_filepaths[file.exists(expected_filepaths)]
    if (length(existing_files) > 0) {
      rlang::abort(
        c("Output files already exist and `overwrite` is FALSE.",
          "i" = "The following files would be overwritten:",
          "*" = paste(basename(existing_files), collapse = ", "),
          "i" = "To proceed, set `overwrite = TRUE` or remove these files manually."))
    }
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --- 3. Input Gathering and Static Index Separation ---
  input_rasters <- list(tmin = tmin, tmax = tmax, prcp = prcp, tavg = tavg, 
                        srad = srad, mois = mois)
  input_rasters <- purrr::discard(input_rasters, is.null)

  valid_static_indices <- c("warmest_unit", "coldest_unit", "wettest_unit", 
    "driest_unit", "high_rad_unit", "low_rad_unit", "high_mois_unit", "low_mois_unit", 
    "warmest_period", "coldest_period", "wettest_period", "driest_period", 
    "high_mois_period", "low_mois_period")

  static_index_rasters <- dot_args[names(dot_args) %in% valid_static_indices]
  other_args <- dot_args[!names(dot_args) %in% valid_static_indices]

  is_spat_rast <- sapply(static_index_rasters, inherits, "SpatRaster")
  if (length(static_index_rasters) > 0 && !all(is_spat_rast)) {
    rlang::abort("All static indices provided via '...' must be SpatRaster objects.")
  }

  all_input_rasters <- c(input_rasters, static_index_rasters)
  if (length(all_input_rasters) == 0) rlang::abort("No input SpatRaster objects were provided.")

  full_raster_stack <- do.call(c, all_input_rasters)

  # --- 4. The Multi-Stage Decision Logic ---
  use_terra_workflow <- FALSE
  if (method == "terra") {
    use_terra_workflow <- TRUE
    rlang::inform("User forced 'terra' workflow.")
  } else if (method == "tiled") {
    use_terra_workflow <- FALSE
    rlang::inform("User forced 'tiled' workflow.")
  } else { # "auto"
    rlang::inform("Using 'auto' method to select workflow...")
    if (is.null(user_region)) {
      check_terra_workflow <- terra::mem_info(full_raster_stack[[1]], 
                                              n = length(full_raster_stack) * terra::nlyr(full_raster_stack[[1]]),
                                              print = FALSE)
      use_terra_workflow <- check_terra_workflow["fits_mem"] == 1
      if (use_terra_workflow) {
        rlang::inform("Full rasters appear to fit in memory. Selecting 'terra' workflow.")
      } else {
        rlang::inform("Full rasters are too large for memory. Selecting 'tiled' workflow.")
      }
    } else {
      template_rast <- full_raster_stack[[1]][[1]]
      if (!terra::relate(terra::ext(user_region), terra::ext(template_rast), "intersects")) {
        rlang::abort("The provided 'user_region' does not overlap with the input rasters.")
      }
      prop_raster <- terra::crop(template_rast, terra::ext(user_region))
      mem_needed_gb <- terra::mem_info(prop_raster, 
                                       n = length(full_raster_stack) * terra::nlyr(full_raster_stack[[1]]),
                                       print = FALSE)
      if (mem_needed_gb["fits_mem"] == 1) {
        use_terra_workflow <- TRUE
        rlang::inform("Estimated cropped region appears to fit in memory. Selecting 'terra' workflow.")
      } else {
        use_terra_workflow <- FALSE
        rlang::inform("Estimated cropped region is likely too large for memory. Selecting 'tiled' workflow.")
      }
    }
  }

  # --- 5. Warn about Misused Arguments ---
  if (use_terra_workflow) {
    if (!missing(tile_degrees)) {
      rlang::warn("The 'tile_degrees' argument is ignored because the 'terra' workflow was selected.")
    }
  }

  # --- 6. Delegate to the Correct Workflow ---
  if (use_terra_workflow) {
    final_rasters <- all_input_rasters
    if (!is.null(user_region)) {
      rlang::inform("Performing crop operation for 'terra' workflow...")
      final_rasters <- lapply(all_input_rasters, function(r) terra::crop(r, user_region, mask = TRUE))
    }
    call_args_terra <- c(list(bios = bios, output_dir = output_dir, period_length = period_length, circular = circular,
        gdal_opt = gdal_opt, overwrite = overwrite), final_rasters, other_args)
    bioclim_results <- do.call(bioclim_terra, call_args_terra)
  } else { # Tiled workflow
    if (any(sapply(all_input_rasters, terra::inMemory))) {
      rlang::abort(
        c("The 'tiled' workflow requires all input SpatRasters (including static indices) to point to files on disk.",
          "i" = "Please save any in-memory rasters to disk first or use `method = 'terra'` if they are small enough."))
    }
    rlang::inform("Extracting file paths from SpatRasters for tiled workflow.")
    input_paths <- lapply(input_rasters, terra::sources)
    static_index_paths <- list()
    if (length(static_index_rasters) > 0) {
      for (i in seq_along(static_index_rasters)) {
        name <- names(static_index_rasters)[i]
        path <- terra::sources(static_index_rasters[[i]])
        static_index_paths[[paste0(name, "_path")]] <- path
      }
    }
    n_units <- terra::nlyr(input_rasters[[1]])
    call_args_tiled <- c(list(bios = bios, n_units = n_units, period_length = period_length, circular = circular,
                              user_region = user_region, tile_degrees = tile_degrees, output_dir = output_dir),
                              input_paths, static_index_paths, other_args)
    temp_results_dir <- do.call(bioclim_fast, call_args_tiled)
    bioclim_results <- write_layers(
      input_dir = temp_results_dir, 
      output_dir = output_dir, 
      file_pattern = "bio",
      gdal_opt = gdal_opt, 
      overwrite = overwrite)
  }
  rlang::inform(paste("Processing complete. Final rasters are in:", normalizePath(output_dir)))
  return(terra::rast(bioclim_results))
}