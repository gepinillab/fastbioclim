#' Derive Custom Summary Statistics from Climate Variables
#'
#' Calculates a wide range of custom summary statistics for a primary climate
#' variable, with options for interactions with a second variable. This function
#' serves as a smart wrapper that automatically selects the most efficient
#' processing workflow (in-memory vs. tiled).
#'
#' @details
#' This function provides a flexible alternative to `derive_bioclim()` for any
#' multi-layer climate variable (e.g., wind speed, humidity). It unifies two
#' processing backends, controlled by the `method` argument:
#' \itemize{
#'   \item `"auto"`: (Default) Intelligently chooses between "terra" and "tiled".
#'   \item `"terra"`: Forces the fast, in-memory workflow.
#'   \item `"tiled"`: Forces the memory-safe, out-of-core workflow. Requires
#'     that all input SpatRasters point to files on disk.
#' }
#'
#' @section Static Indices:
#' For advanced control, provide pre-calculated index rasters as named
#' `SpatRaster` objects via the `...` argument (e.g., `max_unit = max_idx_rast`).
#' Supported indices: `max_unit`, `min_unit`, `max_period`, `min_period`,
#' `max_interactive`, `min_interactive`.
#'
#' @param variable A `terra::SpatRaster` object for the primary variable.
#' @param stats A character vector of statistics to compute for the primary
#'   variable. Supported: `"mean"`, `"max"`, `"min"`, `"sum"`, `"stdev"`, `"cv_cli"`,
#'   `"max_period"`, `"min_period"`.
#' @param inter_variable (Optional) A `terra::SpatRaster` for an interactive variable.
#' @param inter_stats (Optional) A character vector of interactive statistics to compute.
#'   Requires `inter_variable`. Supported: `"max_inter"`, `"min_inter"`.
#' @param prefix_variable A character string used as the prefix for all output file
#'   names (e.g., `prefix_variable = "wind"` results in "wind_mean.tif", "wind_max.tif").
#' @param suffix_inter_max Character. Suffix for the "max_inter" statistic name. Default: "inter_high".
#' @param suffix_inter_min Character. Suffix for the "min_inter" statistic name. Default: "inter_low".
#' @param output_dir The directory where the final summary rasters will be saved.
#' @param period_length Integer. The number of temporal units defining a "period". Default: 3.
#' @param period_stats Character. The statistic ("mean" or "sum") to summarize
#'   data over each period. Default: "mean".
#' @param circular Logical. If `TRUE` (the default), period calculations wrap around.
#' @param user_region (Optional) An `sf` or `terra::SpatVector` object defining the
#'   area of interest.
#' @param method The processing method. See Details for more information.
#' @param tile_degrees (Tiled method only) The approximate size of processing tiles.
#' @param gdal_opt (Optional) A character vector of GDAL creation options for the
#'   output GeoTIFF files.
#' @param overwrite (Optional) Logical. If `FALSE` (the default), the function will
#'   stop if output files already exist.
#' @param ... Additional arguments, primarily for passing static index `SpatRaster`
#'   objects. See the "Static Indices" section.
#'
#' @return A `terra::SpatRaster` object pointing to the newly created summary rasters.
#' @export
derive_statistics <- function(variable,
  stats = c("mean", "max", "min"),
  inter_variable = NULL,
  inter_stats = NULL,
  prefix_variable = "var",
  output_dir = tempdir(),
  period_length = 3,
  period_stats = "mean",
  circular = TRUE,
  user_region = NULL,
  method = c("auto", "tiled", "terra"),
  tile_degrees = 5,
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE,
  ...) {

  # --- 1. Argument Capture and Initial Setup ---
  method <- match.arg(method)
  dot_args <- list(...)

  # --- 2. Fail-Fast Pre-Check for Existing Files ---
  if (!overwrite) {
    # Generate expected output names based on requested stats
    stat_suffixes <- unique(c(stats))
    if ("max_inter" %in% inter_stats) {
      stat_suffixes <- c(stat_suffixes, suffix_inter_max)
    }
    if ("min_inter" %in% inter_stats) {
      stat_suffixes <- c(stat_suffixes, suffix_inter_min)
    }
    expected_filenames <- paste0(prefix_variable, "_", stat_suffixes, ".tif")
    expected_filepaths <- file.path(output_dir, expected_filenames)
    existing_files <- expected_filepaths[file.exists(expected_filepaths)]
    if (length(existing_files) > 0) {
      rlang::abort(c("Potential output files already exist and `overwrite` is FALSE.",
                     "i" = "Please set `overwrite = TRUE` or check the output directory."))
    }
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --- 3. Input Gathering and Static Index Separation ---
  if(is.null(variable)) rlang::abort("The 'variable' SpatRaster must be provided.")
  input_rasters <- list(variable = variable, inter_variable = inter_variable)
  input_rasters <- purrr::discard(input_rasters, is.null)

  valid_static_indices <- c("max_unit", "min_unit", "max_period", "min_period",
                            "max_interactive", "min_interactive")
  static_index_rasters <- dot_args[names(dot_args) %in% valid_static_indices]
  other_args <- dot_args[!names(dot_args) %in% valid_static_indices]
  
  is_spat_rast <- sapply(static_index_rasters, inherits, "SpatRaster")
  if (length(static_index_rasters) > 0 && !all(is_spat_rast)) {
    rlang::abort("All static indices provided via '...' must be SpatRaster objects.")
  }

  all_input_rasters <- c(input_rasters, static_index_rasters)
  full_raster_stack <- do.call(c, all_input_rasters)

  # --- 4. "Auto" Decision Logic ---
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
      final_rasters <- lapply(all_input_rasters, function(r) terra::crop(r, user_region, mask = TRUE))
    }

    # The wrapper handles writing the files
    # Prepare arguments for the internal terra function
    call_args_terra <- c(
      list(stats = stats, inter_stats = inter_stats, prefix_variable = prefix_variable,
           period_length = period_length, period_stats = period_stats, circular = circular,
           output_dir = output_dir, overwrite = overwrite, gdal_opt = gdal_opt),
      final_rasters, other_args)

    # Call the internal function to get the results in memory
    stats_rast_stack <- do.call(stats_terra, call_args_terra)
    bioclim_results <- terra::rast(stats_rast_stack)

  } else { # Tiled workflow
    if (any(sapply(all_input_rasters, terra::inMemory))) {
      rlang::abort("Tiled workflow requires all input SpatRasters to be file-backed.")
    }

    # Prepare arguments for the internal tiled function
    input_paths <- lapply(input_rasters, terra::sources)
    static_index_paths <- list()
    if (length(static_index_rasters) > 0) {
      for(i in seq_along(static_index_rasters)) {
        name <- names(static_index_rasters)[i]
        static_index_paths[[paste0(name, "_path")]] <- terra::sources(static_index_rasters[[i]])
      }
    }
    n_units <- terra::nlyr(variable)

    call_args_tiled <- c(
      list(n_units = n_units, stats = stats, inter_stats = inter_stats,
           prefix_variable = prefix_variable, period_length = period_length, period_stats = period_stats,
           circular = circular, user_region = user_region, tile_degrees = tile_degrees, output_dir = output_dir),
      input_paths, static_index_paths, other_args
    )

    # Create a temporary directory inside the final output directory
    stats_dir <- file.path(output_dir, paste0("stats_qs_", basename(tempfile(pattern = ""))))
    call_args_tiled$stats_dir <- stats_dir
    # Run the tiled process
    stats_results_dir <- do.call(stats_fast, call_args_tiled)
    # Assemble the results using write_layers, with the crucial file_pattern argument
    created_files <- write_layers(
      input_dir = stats_results_dir, 
      output_dir = output_dir,
      file_pattern = prefix_variable,
      gdal_opt = gdal_opt, 
      overwrite = overwrite
    )
    bioclim_results <- terra::rast(created_files)
  }
  rlang::inform(paste("Processing complete. Final rasters are in:", normalizePath(output_dir)))
  return(bioclim_results)
}