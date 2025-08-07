# HELPER FUNCTIONS
# -----------------------

# TERRA mode

# Average Temperature
#' @keywords internal
t_avg <- function(tmin, tmax) {
  misqua <- (tmin + tmax) / 2
  names(misqua) <- paste0("tavg_", 1:terra::nlyr(misqua))
  return(misqua)
}

# CV terra
#' @keywords internal
cv_cli <- function(prcp) {
  pr <- prcp + 1
  x <- terra::mean(abs(pr))
  cv <- 100 * terra::stdev(pr, pop = FALSE, na.rm = TRUE) / x
  return(cv)
}

#' Windows
#'
#' @param x spatRaster
#' @param period Length of period. Default is three. If you are using months. It will be a quarter.
#' @param circular logical Include first month/weeks?
#' @keywords internal
get_window <- function(x, period, circular)  {
  lng <- terra::nlyr(x)
  if (circular == TRUE) {
    ind <- c(1:lng,  1:(period - 1))
    m <- matrix(ncol = period, nrow = lng)
    for (i in 1:period) {
      m[, i] <- ind[i:(lng + i - 1)]
      # if (i != 3) {
      #   m[, i] <- ind[i:(lng + i - 1)]
      # } else {
      #   m[, i] <- ind[(i - 1):(lng + i - 2)]
      # }
    }
  }
  if (circular == FALSE) {
    ind <- c(1:lng)
    m <- matrix(ncol = period, nrow = lng - period + 1)
    for (i in 1:period) {
      m[, i] <- ind[i:(lng - period + i)]
    }
  }
  vent <- NULL
  for (j in 1:nrow(m)) {
    sum_period <- terra::app(x[[m[j, ]]], sum, na.rm = TRUE)
    vent <- c(vent, sum_period)
  }
  return(terra::rast(vent))
}

#' @keywords internal
testGeom <- function(x, y) {
  testGeom <- terra::compareGeom(x, y, lyrs = TRUE)
  if (testGeom == TRUE) {
    return(x)
  }
}

#' @keywords internal
mismatch_NA <- function(layer) {
  # Get number of layers
  num_lyr <- terra::nlyr(layer)
  # Create a raster that sums if cells has values for all of them
  sum_lyr <- sum(!is.na(layer))
  # Create vector to check if there is a mismatch between NAs
  v_unique <- unique(as.vector(sum_lyr))
  # Check for cell values lower than the number of layers. This means that there
  # are cells with NA values
  if (any(v_unique != 0 & v_unique != num_lyr)) {
    miss_na <- v_unique[v_unique != 0 & v_unique != num_lyr]
    cells_na <- unlist(terra::cells(sum_lyr, miss_na))
    layer_na <- sapply(as.list(cells_na),
                       function(x) {which(is.na(extract(layer, x)))})
    message(paste0("Unexpected NA value in '", substitute(layer), "' object",
                   " | Layer number: ", layer_na,
                   " | Cell ID: ", cells_na, "\n"))

  }
  return(list(logical = any(v_unique != 0 & v_unique != num_lyr),
              sum_lyr = sum_lyr))
}

#' @title Print Bioclimatic Variable Names
#' @description
#' This function prints the names of bioclimatic variables based on the specified indices.
#' @param bios Numeric vector indicating the indices of bioclimatic variables to print.
#' Default is 1:35, which prints all variable names.
#' @return None. Prints the names of the selected bioclimatic variables to the console.
#' @examples
#' bionames()           # Print all bioclimatic variable names
#' bionames(c(1, 5, 12)) # Print names for variables 1, 5, and 12
#' @export
bionames <- function(bios = 1:35) {
  # Bioclimatic variable names
  bioclim_vars <- c(
    "bio01: Annual Mean Temperature",
    "bio02: Mean Diurnal Range",
    "bio03: Isothermality",
    "bio04: Temperature Seasonality",
    "bio05: Max Temperature of Warmest Period",
    "bio06: Min Temperature of Coldest Period",
    "bio07: Temperature Annual Range",
    "bio08: Mean Temperature of Wettest Period",
    "bio09: Mean Temperature of Driest Period",
    "bio10: Mean Temperature of Warmest Period",
    "bio11: Mean Temperature of Coldest Period",
    "bio12: Annual Precipitation",
    "bio13: Precipitation of Wettest Period",
    "bio14: Precipitation of Driest Period",
    "bio15: Precipitation Seasonality",
    "bio16: Precipitation of Wettest Period",
    "bio17: Precipitation of Driest Period",
    "bio18: Precipitation of Warmest Period",
    "bio19: Precipitation of Coldest Period",
    "bio20: Annual Mean Radiation",
    "bio21: Highest Period Radiation",
    "bio22: Lowest Period Radiation",
    "bio23: Radiation Seasonality",
    "bio24: Radiation of Wettest Period",
    "bio25: Radiation of Driest Period",
    "bio26: Radiation of Warmest Period",
    "bio27: Radiation of Coldest Period",
    "bio28: Annual Mean Moisture Content",
    "bio29: Highest Period Moisture Content",
    "bio30: Lowest Period Moisture Content",
    "bio31: Moisture Content Seasonality",
    "bio32: Mean Moisture Content of Most Moist Period",
    "bio33: Mean Moisture Content of Least Moist Period",
    "bio34: Mean Moisture Content of Warmest Period",
    "bio35: Mean Moisture Content of Coldest Period"
  )
  
  # Validate input
  if (!all(bios %in% 1:35)) {
    stop("The 'bios' parameter must contain only numbers between 1 and 35.")
  }
  
  # Print the selected variable names
  cat(bioclim_vars[bios], sep = "\n")
}

# ======================================================================

# FAST mode

# --- Core Helper Functions ---

`%||%` <- function(x, y) { if (is.null(x)) y else x }

#' Standardize Environmental Variable Input Format
#'
#' Converts data frames, matrices, or vectors into a consistent matrix format
#' with preserved column/row names for downstream processing.
#'
#' @param evar Input environmental variable data. Can be:
#'   - data.frame: Converted to matrix
#'   - vector: Transposed to row matrix
#'   - matrix: Preserved with names checked
#'
#' @return Matrix with column names preserved. If input was a vector,
#'         returns 1-row matrix with vector names as column names.
#' @keywords internal
check_evar <- function(evar){
  names_col <- colnames(evar)
  if(is.data.frame(evar)){
    evar <- Rfast::data.frame.to_matrix(evar)
  } else if(is.vector(evar)){
    names_col <- names(evar)
    evar <- t(evar) # Transpose vector to a 1-row matrix
  }
  # Add check if it's already a matrix? - Assumed handled if not df or vector
  colnames(evar) <- names_col

  return(evar)
}

#' Calculate Sliding Periods for Temporal Analysis
#'
#' Defines sliding periods of a specified length for temporal calculations,
#' handling wrap-around cases for circular data.
#'
#' @param n_units Integer. The total number of temporal units (layers).
#' @param period_length Integer. The number of consecutive units in each period.
#' @param circular Logical. If TRUE, allows periods to wrap around.
#' @return List where each element contains `period_length` consecutive unit indices.
#' @keywords internal
compute_periods <- function(n_units, period_length, circular = TRUE) {
  if (period_length <= 0 || period_length > n_units) stop("'period_length' invalid.")
  if (n_units <= 0) stop("'n_units' must be positive.")

  indices <- seq_len(n_units)
  num_periods <- ifelse(circular, n_units, n_units - period_length + 1)
  periodos <- vector("list", num_periods)

  for (i in seq_len(num_periods)) {
    start_index <- i
    current_period_indices <- numeric(period_length)
    for (j in seq_len(period_length)) {
      current_index <- start_index + j - 1
      if (current_index > n_units && circular) {
        current_period_indices[j] <- indices[((current_index - 1) %% n_units) + 1]
      } else if (current_index <= n_units) {
        current_period_indices[j] <- indices[current_index]
      } else { stop("Indexing error in compute_periods (non-circular).") }
    }
     periodos[[i]] <- current_period_indices
  }
  return(periodos)
}


#' Calculate Temporal Period Aggregates
#'
#' Computes aggregates (sums or mean) over defined temporal periods and identifies min/max periods.
#'
#' @param variable Matrix of temporal values (rows=pixels/cells, cols=units).
#' @param periodos List of unit groupings (output from compute_periods).
#' @param n_units Integer. Total number of temporal units.
#' @param period_length Integer. Length of each period.
#' @param stat Character. Either `"mean"` or `"sum"`, specifying whether to calculate the average or total
#' @return Matrix with period sums (`num_periods` columns named "period_X"),
#'   "min_idx", and "max_idx" columns.
#' @keywords internal
var_periods <- function(variable, periodos, n_units, period_length, stat) {
  if (length(stat) != 1 || !(stat %in% c("mean", "sum"))) {
    stop('`stat` must be exactly one of "mean" or "sum".')
  }
  
  num_periods_calculated <- length(periodos)
  
  pnames_base <- paste0("period_", seq_len(num_periods_calculated))
  
  pnames <- c(pnames_base, "min_idx", "max_idx")
  
  vperiods <- sapply(periodos, function(p_indices) {
    varval <- variable[, p_indices, drop = FALSE]
    if (stat == "mean") {
      return(Rfast::rowmeans(varval))
    } else {
      return(Rfast::rowsums(varval))
    }
  })
  
    # Handle case where input is a single row
    if (is.vector(vperiods)) vperiods <- t(vperiods)
  
    # Get indices (not values) of min and max period per row
    pmin <- Rfast::rowMins(vperiods, value = FALSE)
    pmax <- Rfast::rowMaxs(vperiods, value = FALSE)
    
    vperiods <- cbind(vperiods, pmin, pmax)
    colnames(vperiods) <- pnames
    
    return(vperiods)
}

#' Create an Inverse Cell ID Translation Function
#'
#' Generates a function to translate cell IDs from a source raster grid
#' to a target raster grid, considering potential offsets and different dimensions.
#' Handles cases where the target grid is a subset (e.g., cropped/masked)
#' of the source grid.
#'
#' @param ncol_src Integer. Number of columns in the source raster.
#' @param ncol_tgt Integer. Number of columns in the target raster.
#' @param row_offset Integer. Row offset of the target grid's top-left corner
#'   relative to the source grid's top-left corner (0-based).
#' @param col_offset Integer. Column offset of the target grid's top-left corner
#'   relative to the source grid's top-left corner (0-based).
#'
#' @return A function that takes a vector of source cell IDs (`cell_src`) and
#'   returns a vector of corresponding target cell IDs. Cells falling outside
#'   the target grid bounds will have `NA_integer_` as their target ID.
#' @keywords internal
define_translate <- function(ncol_src, ncol_tgt, row_offset, col_offset) {
  force(ncol_src) # Force evaluation of arguments in the enclosing environment
  force(ncol_tgt)
  force(row_offset)
  force(col_offset)

  function(cell_src) {
    if (!is.integer(cell_src)) cell_src <- as.integer(cell_src) # Ensure integer input

    # Preallocate result vector with NA
    result <- rep(NA_integer_, length(cell_src))

    # Calculate source row and column (1-based)
    row_src <- ((cell_src - 1L) %/% ncol_src) + 1L
    col_src <- ((cell_src - 1L) %% ncol_src) + 1L

    # Calculate potential target row and column (1-based) by applying offset
    row_tgt <- row_src - row_offset
    col_tgt <- col_src - col_offset

    # Identify valid cells (within the bounds of the target grid)
    valid <- row_tgt >= 1L & col_tgt >= 1L

    # Compute target cell ID *only* for valid indices
    result[valid] <- (row_tgt[valid] - 1L) * ncol_tgt + col_tgt[valid]
    result
  }
}
