#' Validate Climate Input Rasters
#'
#' Performs rigorous checks on a set of climate raster files to ensure they
#' are suitable for bioclimatic variable calculation.
#'
#' @details
#' This function checks for two main sources of error:
#' 1.  **Geometric Inconsistency**: Ensures all input rasters share the exact same
#'     coordinate reference system (CRS), extent, and resolution.
#' 2.  **NA Mismatch**: Checks if the pattern of NA values is consistent across
#'     all layers of all provided climate variables. Mismatched NAs can lead
#'     to silent errors in calculations. This check can be time-consuming for
#'     very large rasters.
#'
#' @param ... Named arguments providing character vectors of file paths for each
#'   climate variable (e.g., `tmin_files = c(...)`, `prcp_files = c(...)`).
#' @param check_nas Logical. If `TRUE` (the default), perform the potentially
#'   slow check for mismatched NA values. If `FALSE`, only check for geometry.
#' @return A list containing a summary of the validation:
#'   \item{`is_valid`}{A single logical value: `TRUE` if all checks pass, `FALSE` otherwise.}
#'   \item{`geom_report`}{A message detailing the result of the geometry check.}
#'   \item{`na_report`}{A message detailing the result of the NA check.}
#'   \item{`mismatch_raster`}{If `check_nas=TRUE` and mismatches are found, a `SpatRaster` where non-zero values indicate pixels with inconsistent NAs.}
#' @export
#' @importFrom terra rast compareGeom tapp lapp
#' @importFrom purrr discard reduce
check_rasters <- function(..., check_nas = TRUE) {
  
  files_list <- list(...)
  if (length(files_list) == 0) {
    stop("No input files provided.")
  }

  message("Creating SpatRaster objects for validation...")
  raster_list <- lapply(files_list, terra::rast)
  
  report <- list(
    is_valid = TRUE,
    geom_report = "Geometry check not performed.",
    na_report = "NA check not performed.",
    mismatch_raster = NULL
  )

  # 1. Geometry Check (fast)
  tryCatch({
    purrr::reduce(raster_list, terra::compareGeom)
    report$geom_report <- "SUCCESS: All rasters share the same geometry (CRS, extent, resolution)."
  }, error = function(e) {
    report$is_valid <<- FALSE
    report$geom_report <<- paste("FAILURE: Geometries are inconsistent.", e$message)
    # Stop now, no point checking NAs if geometry differs
    return(report)
  })

  if (!report$is_valid) return(report)
  
  # 2. NA Value Check (can be slow)
  if (check_nas) {
    message("Performing NA consistency check (this may take a while)...")
    
    # For each stack, create a single layer summary: 0=all NA, 1=some NA, 2=no NA
    na_summaries <- lapply(raster_list, function(r) {
      terra::lapp(r, fun = function(x) {
        if (all(is.na(x))) return(0)
        if (any(is.na(x))) return(1)
        return(2)
      })
    })
    
    # Stack the summaries and check for pixel-wise inconsistencies
    summary_stack <- terra::rast(na_summaries)
    
    # A pixel is inconsistent if the summary values (0,1,2) are not all the same
    mismatch_rast <- terra::app(summary_stack, fun = function(x) {
      # The range (max-min) will be > 0 if values are not identical
      max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
    })
    
    mismatch_count <- terra::global(mismatch_rast, "sum", na.rm = TRUE)$sum
    
    if (mismatch_count > 0) {
      report$is_valid <- FALSE
      report$na_report <- sprintf("FAILURE: Found %d pixels with inconsistent NA patterns across variables.", mismatch_count)
      report$mismatch_raster <- mismatch_rast
    } else {
      report$na_report <- "SUCCESS: NA patterns are consistent across all variables."
    }
  }

  message("Validation complete.")
  return(report)
}