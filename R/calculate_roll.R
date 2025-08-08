#' Calculate Rolling Temporal Averages for SpatRasters
#'
#' Calculates temporal summaries for each time unit over a moving window of cycles.
#' This function is designed for time series where fundamental time **units** (e.g., months)
#' are grouped into repeating **cycles** (e.g., years).
#'
#' @param x A `terra::SpatRaster` object where each layer represents a time **unit**.
#' @param window_size Integer. The size of the moving window, measured in the number of **cycles**.
#'   For example, if the data cycle is annual (`freq = 12`), a `window_size` of 20
#'   represents a 20-year window.
#' @param freq Integer. The number of time **units** (layers) that constitute one
#'   complete **cycle**. Common examples: 12 for monthly units in a yearly cycle,
#'   or 24 for hourly units in a daily cycle.
#' @param step Integer. The number of **cycles** to slide the window by for each
#'   iteration. Default is 1.
#' @param fun Character. The name of the summary function (e.g., "mean"). Default is "mean".
#' @param output_prefix A character string for output filenames. Default is "output".
#' @param output_dir Directory to save the final GeoTIFF files.
#' @param name_template A character string defining the template for output filenames,
#'   using `glue` syntax. Default: `"{prefix}_w{start_window}-{end_window}_u{idx_unit}"`.
#'   Available placeholders are:
#'   \itemize{
#'     \item `{prefix}`: The value from `output_prefix`.
#'     \item `{start_window}`: The starting **cycle** index of the window.
#'     \item `{end_window}`: The ending **cycle** index of the window.
#'     \item `{idx_unit}`: The index of the time **unit** within the cycle (e.g., the month number).
#'   }
#' @param user_region (Optional) An `sf` or `terra::SpatVector` for clipping.
#' @param method Processing method: "auto", "tiled", or "terra".
#' @param tile_degrees (Tiled method only) Approximate size of processing tiles.
#' @param gdal_opt (Optional) GDAL creation options for GeoTIFFs.
#' @param overwrite Logical. If `FALSE` (default), stops if output files exist.
#'
#' @return A `terra::SpatRaster` object pointing to the newly created files.
#' @export
calculate_roll <- function(x,
  window_size,
  freq = 12,
  step = 1,
  fun = "mean",
  name_template = "{prefix}_w{start_window}-{end_window}_u{idx_unit}",
  output_prefix = "output",
  output_dir = tempdir(),
  user_region = NULL,
  method = c("auto", "tiled", "terra"),
  tile_degrees = 5,
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE) {

  # --- 1. Argument Validation and Setup ---
  method <- match.arg(method)
  if (!inherits(x, "SpatRaster")) rlang::abort("Input 'x' must be a SpatRaster.")
  n_in <- nlyr(x)

  # Validate time series parameters
  if (missing(window_size)) rlang::abort("'window_size' is required.")
  if (n_in %% freq != 0) {
    rlang::abort(glue::glue("Total layers ({n_in}) is not divisible by frequency ({freq})."))
  }
  total_cycles <- n_in / freq
  if (window_size > total_cycles) {
    rlang::abort(glue::glue("window_size ({window_size}) cannot be larger than total cycles ({total_cycles})."))
  }

  start_units <- seq(1, total_cycles - window_size + 1, by = step)

  output_names_list <- lapply(start_units, function(start_w) {
    end_w <- start_w + window_size - 1
    
    sapply(1:freq, function(p) {
      fun <- fun
      prefix <- output_prefix
      start_window <- sprintf(paste0("%0", nchar(end_w), "d"), start_w)
      end_window <- end_w
      idx_unit <- sprintf(paste0("%0", nchar(freq), "d"), p)
      
      # Aplicamos la plantilla
      glue::glue(name_template)
    })
  })
  
  all_output_names <- unlist(output_names_list)

  # Flatten the list to get all expected output names
  all_output_names <- unlist(output_names_list)
  n_units_out <- length(all_output_names)
  file_pattern_prefix <- if (is.null(output_prefix)) fun else output_prefix

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --- 2. Fail-Fast Pre-Check for Existing Files ---
  if (!overwrite) {
    expected_filepaths <- file.path(output_dir, paste0(all_output_names, ".tif"))
    if (any(file.exists(expected_filepaths))) {
      rlang::abort(
      c("Output files already exist and `overwrite` is FALSE.",
      "i" = "To proceed, set `overwrite = TRUE` or remove the existing files.")
      )
    }
  }

  # --- 3. The "Auto" Decision Logic ---
  use_terra_workflow <- FALSE
  if (method == "terra") {
    use_terra_workflow <- TRUE
    rlang::inform("User forced 'terra' (in-memory) workflow.")
  } else if (method == "tiled") {
    use_terra_workflow <- FALSE
    rlang::inform("User forced 'tiled' (out-of-core) workflow.")
  } else { # "auto"
    rlang::inform("Using 'auto' method to select workflow...")
    # For rolling average, the memory peak is one window + its results
    n_layers_for_check <- (window_size * freq) + freq
    template_rast <- if (is.null(user_region)) x[[1]] else crop(x[[1]], user_region)
    mem_info <- terra::mem_info(template_rast, n = n_layers_for_check, print = FALSE)
    if (mem_info["fits_mem"] == 1) {
      use_terra_workflow <- TRUE
      rlang::inform("A single window appears to fit in memory. Selecting 'terra' workflow.")
    } else {
      use_terra_workflow <- FALSE
      rlang::inform("A single window may be too large for memory. Selecting 'tiled' workflow.")
    }
  }

  # --- 4. Delegate to the Correct Workflow ---
  if (use_terra_workflow) {
    final_raster <- if (!is.null(user_region)) {
    rlang::inform("Performing crop operation for 'terra' workflow...")
    terra::crop(x, user_region, mask = TRUE)
    } else {
    x
    }

    result_paths <- roll_terra(
      x = final_raster,
      window_size = window_size,
      freq = freq,
      step = step,
      fun = fun,
      output_names_list = output_names_list, # Pass the list of names
      output_dir = output_dir,
      gdal_opt = gdal_opt,
      overwrite = overwrite
    )

  } else { # Tiled workflow
    if (all(terra::inMemory(x))) {
      rlang::abort("The 'tiled' workflow requires the input SpatRaster to point to a file on disk.")
    }

    temp_qs_dir <- roll_fast(
      paths = terra::sources(x),
      window_size = window_size,
      freq = freq,
      step = step,
      fun = fun,
      output_names_list = output_names_list,
      user_region = user_region,
      tile_degrees = tile_degrees,
      output_dir = output_dir
    )

    # Use write_layers to assemble the results
    result_paths <- write_layers(
      input_dir = temp_qs_dir,
      output_dir = output_dir,
      file_pattern = file_pattern_prefix,
      gdal_opt = gdal_opt,
      overwrite = overwrite
    )
  }
  rlang::inform(paste("Processing complete. Final rasters are in:", normalizePath(output_dir)))
  return(terra::rast(result_paths))
}