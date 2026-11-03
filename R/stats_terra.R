#' In-Memory Custom Variable Summarization
#'
#' Internal function to calculate custom summary statistics using `terra` functions.
#' Designed for datasets that fit into RAM.
#'
#' @param variable A SpatRaster for the primary variable.
#' @param stats Character vector of stats to compute.
#' @param period_length Integer, length of a calculation period.
#' @param period_stats Character, either `"mean"` or `"sum"`; statistic used to
#'   summarize each period.
#' @param circular Logical, whether to wrap periods.
#' @param inter_variable (Optional) A `terra::SpatRaster` for an interactive
#'   variable.
#' @param inter_stats Character vector of interactive statistics to compute
#'   (`"max_inter"`, `"min_inter"`).
#' @param max_unit,min_unit,max_period,min_period,max_interactive,min_interactive
#'   (Optional) `terra::SpatRaster` with a static index.
#' @param output_prefix Character, prefix for output layer names.
#' @param suffix_inter_max,suffix_inter_min Character, suffixes for the
#'   `"max_inter"` and `"min_inter"` statistic names.
#' @param gdal_opt Character vector of GDAL options for writing.
#' @param overwrite Logical, whether to overwrite existing files.
#' @param output_dir Character, path to save final rasters.
#' @param verbose Logical, If `TRUE`, prints messages.
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
                        output_prefix = "var",
                        suffix_inter_max = "inter_high", 
                        suffix_inter_min = "inter_low",
                        gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
                        overwrite = FALSE, 
                        output_dir = tempdir(),
                        verbose = TRUE,
                        ...) {
  # Check for same extent, number of rows and columns, projection,
  # resolution, and origin
  sameGeom <- class(purrr::reduce(list(variable, inter_variable, max_unit, min_unit, 
                                       max_period, min_period, 
                                       max_interactive, min_interactive) |>
                                    purrr::discard(is.null),
                                  testGeom))
  if (sameGeom == "SpatRaster") {
    if (verbose) message("SpatRasters have same extent, number of rows and columns, ",
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
    names(mean_stat) <- paste0(output_prefix, "_mean")
  }
  
  # SUM
  if ("sum" %in% stats) {
    sum_stat <- terra::app(variable, "sum", na.rm = TRUE)
    names(sum_stat) <- paste0(output_prefix, "_sum")
  }
  
  # MAX
  if ("max" %in% stats & is.null(max_unit)) {
    max_stat <- terra::app(variable, "max", na.rm = TRUE)
    names(max_stat) <- paste0(output_prefix, "_max")
  } else if ("max" %in% stats & !is.null(max_unit)) {
    max_stat <- terra::selectRange(variable, max_unit)
    names(max_stat) <- paste0(output_prefix, "_max")
  }
  
  # MIN
  if ("min" %in% stats & is.null(min_unit)) {
    min_stat <- terra::app(variable, "min", na.rm = TRUE)
    names(min_stat) <- paste0(output_prefix, "_min")
  } else if ("min" %in% stats & !is.null(min_unit)) {
    min_stat <- terra::selectRange(variable, min_unit)
    names(min_stat) <- paste0(output_prefix, "_min")
  }
  
  # STDEV
  if ("stdev" %in% stats) {
    # Sample standard deviation, as in bio04/bio31 and stats_fast()
    stdev_stat <- terra::stdev(variable, pop = FALSE, na.rm = TRUE)
    names(stdev_stat) <- paste0(output_prefix, "_stdev")
  }
  
  # CV_CLI
  if ("cv_cli" %in% stats) {
    cv_cli_stat <- cv_cli(variable)
    names(cv_cli_stat) <- paste0(output_prefix, "_cv")
  }
  
  # MAX & MIN PERIOD
  if (any(c("max_period", "min_period", "max_inter", "min_inter") %in% c(stats, inter_stats))) {
    # Get windows
    period_windows <- get_window(variable, period_length, circular)
    if (period_stats == "mean") period_windows <- period_windows / period_length
  }
  ## MAX_PERIOD
  if ("max_period" %in% stats & is.null(max_period)) {
    max_period_stat <- terra::app(period_windows, "max", na.rm = TRUE)
    names(max_period_stat) <- paste0(output_prefix, "_max_period")
  } else if ("max_period" %in% stats & !is.null(max_period)) {
    max_period_stat <- terra::selectRange(period_windows, max_period)
    names(max_period_stat) <- paste0(output_prefix, "_max_period")
  }
  ## MIN_PERIOD
  if ("min_period" %in% stats & is.null(min_period)) {
    min_period_stat <- terra::app(period_windows, "min", na.rm = TRUE)
    names(min_period_stat) <- paste0(output_prefix, "_min_period")
  } else if ("min_period" %in% stats & !is.null(min_period)) {
    min_period_stat <- terra::selectRange(period_windows, min_period)
    names(min_period_stat) <- paste0(output_prefix, "_min_period")
  }
  
  # INTERACTIVES (MAX & MIN)
  if (!is.null(inter_variable)) {
    if (("max_inter" %in% inter_stats & is.null(max_interactive)) | 
        ("min_inter" %in% inter_stats & is.null(min_interactive))) {
      inter_windows <- get_window(inter_variable, period_length, circular)
    }
    if ("max_inter" %in% inter_stats) {
      if (is.null(max_interactive)) {
        max_inter_stat <- terra::selectRange(period_windows, terra::which.max(inter_windows))
        names(max_inter_stat) <- paste0(output_prefix, "_", suffix_inter_max)
      } else {
        max_inter_stat <- terra::selectRange(period_windows, max_interactive)
        names(max_inter_stat) <- paste0(output_prefix, "_", suffix_inter_max)
      }
    }
    if ("min_inter" %in% inter_stats) {
      if (is.null(min_interactive)) {
        min_inter_stat <- terra::selectRange(period_windows, terra::which.min(inter_windows))
        names(min_inter_stat) <- paste0(output_prefix, "_", suffix_inter_min)
      } else {
        min_inter_stat <- terra::selectRange(period_windows, min_interactive)
        names(min_inter_stat) <- paste0(output_prefix, "_", suffix_inter_min)
      }
    }
  }
  
  # NA policy: every statistic is NA where any unit of `variable` is NA; the
  # interactive statistics also where any unit of `inter_variable` is NA, unless
  # a static index replaces it. stats_fast() applies the same rule.
  na_var <- anyNA(variable)
  mask_na <- function(r, na) if (is.null(r)) NULL else terra::mask(r, na, maskvalues = TRUE)
  mean_stat <- mask_na(mean_stat, na_var)
  sum_stat <- mask_na(sum_stat, na_var)
  max_stat <- mask_na(max_stat, na_var)
  min_stat <- mask_na(min_stat, na_var)
  stdev_stat <- mask_na(stdev_stat, na_var)
  cv_cli_stat <- mask_na(cv_cli_stat, na_var)
  max_period_stat <- mask_na(max_period_stat, na_var)
  min_period_stat <- mask_na(min_period_stat, na_var)
  na_inter <- if (is.null(inter_variable)) na_var else na_var | anyNA(inter_variable)
  max_inter_stat <- mask_na(max_inter_stat, if (is.null(max_interactive)) na_inter else na_var)
  min_inter_stat <- mask_na(min_inter_stat, if (is.null(min_interactive)) na_inter else na_var)

  raster_stats <- terra::rast(list(
    mean_stat, sum_stat, max_stat, min_stat, stdev_stat, cv_cli_stat, 
    max_period_stat, min_period_stat, max_inter_stat, min_inter_stat) |>
      purrr::discard(is.null))

  
  output_files <- file.path(output_dir, paste0(names(raster_stats), ".tif"))
  
  if (verbose) message("Writing GeoTIFFs...")
  terra::writeRaster(
    raster_stats,
    filename = output_files,
    overwrite = overwrite,
    gdal = gdal_opt
  )
  return(output_files)
}
