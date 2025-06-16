#' In-Memory Custom Variable Summarization
#'
#' Internal function to calculate custom summary statistics using `terra` functions.
#' Designed for datasets that fit into RAM.
#'
#' @param variable A SpatRaster for the primary variable.
#' @param stats Character vector of stats to compute.
#' @param prefix_variable Character, prefix for output layer names.
#' @param ... Other arguments including inter_variable, period_length, circular,
#'   static index SpatRasters, etc.
#' @return A `terra::SpatRaster` object containing the calculated summary layers.
#' @keywords internal
#' @seealso The user-facing wrapper function `derived_statistics()`.
stats_terra <- function(variable, 
                        stats = c("mean", "max", "min", 
                                  "cv_cli", "max_period", "min_period"),
                        period_length = 3,
                        period_stats = "mean",
                        circular = TRUE,
                        inter_variable = NULL,
                        inter_stats = c("max_inter", "min_inter"),
                        max_unit = NULL,
                        min_unit = NULL,
                        max_period = NULL,
                        min_period = NULL,
                        max_interactive = NULL,
                        min_interactive = NULL,
                        prefix_variable = "var",
                        suffix_inter_max = "inter_high", 
                        suffix_inter_min = "inter_low",
                        gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
                        overwrite = FALSE, 
                        output_dir = tempdir(),
                        ...) {
  # Check for same extent, number of rows and columns, projection,
  # resolution, and origin
  sameGeom <- class(purrr::reduce(list(variable, inter_variable, max_unit, min_unit, 
                                       max_period, min_period, 
                                       max_interactive, min_interactive) |>
                                    purrr::discard(is.null),
                                  fastbioclim::testGeom))
  if (sameGeom == "SpatRaster") {
    message("SpatRasters have same extent, number of rows and columns, ",
            "projection, resolution, and origin")
  }
  # Find invalid characters
  inv_stats <- setdiff(stats, c("mean", "max", "min", "sum",
                                "stdev", "cv_cli", "max_period", "min_period"))
  if (length(inv_stats) > 0) {
    stop(paste("Invalid stats character(s) provided:", paste(inv_stats, collapse = ", ")))
  }
  
  if (length(period_stats) != 1 & !is.null(period_stats)) {
    stop(paste("Please provide just one period stats (mean OR sum)"))
  }
  
  inv_period_stats <- setdiff(period_stats, c("mean", "sum"))
  if (length(inv_period_stats) > 0) {
    stop(paste("Invalid period stats character(s) provided:", paste(inv_period_stats, collapse = ", ")))
  }
  
  inv_inter_stats <- setdiff(inter_stats, c("max_inter", "min_inter"))
  if (length(inv_period_stats) > 0) {
    stop(paste("Invalid interactive stats character(s) provided:", paste(inv_inter_stats, collapse = ", ")))
  }
    
  # Create NULL spatRaster
  mean_stat <- sum_stat <- max_stat <- min_stat <- stdev_stat <- cv_cli_stat <- NULL
  max_period_stat <- min_period_stat <- max_inter_stat <- min_inter_stat <- NULL
  
  # MEAN
  if ("mean" %in% stats) {
    mean_stat <- terra::app(variable, "mean", na.rm = TRUE)
    names(mean_stat) <- paste0(prefix_variable, "_mean")
  }
  
  # SUM
  if ("sum" %in% stats) {
    sum_stat <- terra::app(variable, "sum", na.rm = TRUE)
    names(sum_stat) <- paste0(prefix_variable, "_sum")
  }
  
  # MAX
  if ("max" %in% stats & is.null(max_unit)) {
    max_stat <- terra::app(variable, "max", na.rm = TRUE)
    names(max_stat) <- paste0(prefix_variable, "_max")
  } else if ("max" %in% stats & !is.null(max_unit)) {
    max_stat <- terra::selectRange(variable, max_unit)
    names(max_stat) <- paste0(prefix_variable, "_max")
  }
  
  # MIN
  if ("min" %in% stats & is.null(min_unit)) {
    min_stat <- terra::app(variable, "min", na.rm = TRUE)
    names(min_stat) <- paste0(prefix_variable, "_min")
  } else if ("min" %in% stats & !is.null(min_unit)) {
    min_stat <- terra::selectRange(variable, min_unit)
    names(min_stat) <- paste0(prefix_variable, "_min")
  }
  
  # STDEV
  if ("stdev" %in% stats) {
    stdev_stat <- terra::app(variable, "std", na.rm = TRUE)
    names(stdev_stat) <- paste0(prefix_variable, "_stdev")
  }
  
  # CV_CLI
  if ("cv_cli" %in% stats) {
    cv_cli_stat <- fastbioclim::cv_cli(variable)
    names(cv_cli_stat) <- paste0(prefix_variable, "_cv")
  }
  
  # MAX & MIN PERIOD
  if (any(c("max_period", "min_period", "max_inter", "min_inter") %in% c(stats, inter_stats))) {
    # Get windows
    period_windows <- fastbioclim::get_window(variable, period_length, circular)
    if (period_stats == "mean") period_windows <- period_windows / period_length
  }
  ## MAX_PERIOD
  if ("max_period" %in% stats & is.null(max_period)) {
    max_period_stat <- terra::app(period_windows, "max", na.rm = TRUE)
    names(max_period_stat) <- paste0(prefix_variable, "_max_period")
  } else if ("max_period" %in% stats & !is.null(max_period)) {
    max_period_stat <- terra::selectRange(period_windows, max_period)
    names(max_period_stat) <- paste0(prefix_variable, "_max_period")
  }
  ## MIN_PERIOD
  if ("min_period" %in% stats & is.null(min_period)) {
    min_period_stat <- terra::app(period_windows, "min", na.rm = TRUE)
    names(min_period_stat) <- paste0(prefix_variable, "_min_period")
  } else if ("min_period" %in% stats & !is.null(min_period)) {
    min_period_stat <- terra::selectRange(period_windows, min_period)
    names(min_period_stat) <- paste0(prefix_variable, "_min_period")
  }
  
  # INTERACTIVES (MAX & MIN)
  if (!is.null(inter_variable)) {
    if (("max_inter" %in% inter_stats & is.null(max_interactive)) | 
        ("min_inter" %in% inter_stats & is.null(min_interactive))) {
      inter_windows <- fastbioclim::get_window(inter_variable, period_length, circular)
    }
    if ("max_inter" %in% inter_stats) {
      if (is.null(max_interactive)) {
        max_inter_stat <- terra::selectRange(period_windows, terra::which.max(inter_windows))
        names(max_inter_stat) <- paste0(prefix_variable, "_", suffix_inter_max)
      } else {
        max_inter_stat <- terra::selectRange(period_windows, max_interactive)
        names(max_inter_stat) <- paste0(prefix_variable, "_", suffix_inter_max)
      }
    }
    if ("min_inter" %in% inter_stats) {
      if (is.null(min_interactive)) {
        min_inter_stat <- terra::selectRange(period_windows, terra::which.min(inter_windows))
        names(min_inter_stat) <- paste0(prefix_variable, "_", suffix_inter_min)
      } else {
        min_inter_stat <- terra::selectRange(period_windows, min_interactive)
        names(min_inter_stat) <- paste0(prefix_variable, "_", suffix_inter_min)
      }
    }
  }
  
  raster_stats <- terra::rast(list(
    mean_stat, sum_stat, max_stat, min_stat, stdev_stat, cv_cli_stat, 
    max_period_stat, min_period_stat, max_inter_stat, min_inter_stat) |>
      purrr::discard(is.null))

  
  output_files <- file.path(output_dir, paste0(names(raster_stats), ".tif"))
  
  message("Writing GeoTIFFs...")
  terra::writeRaster(
    raster_stats,
    filename = output_files,
    overwrite = overwrite,
    gdal = gdal_opt
  )
  return(output_files)
}
