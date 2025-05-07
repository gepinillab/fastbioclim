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
#' Computes aggregates (sums) over defined temporal periods and identifies min/max periods.
#'
#' @param variable Matrix of temporal values (rows=pixels/cells, cols=units).
#' @param periodos List of unit groupings (output from compute_periods).
#' @param n_units Integer. Total number of temporal units.
#' @param period_length Integer. Length of each period.
#' @return Matrix with period sums (`num_periods` columns named "var_period_X"),
#'   "min_periodID", and "max_periodID" columns.
#' @keywords internal
var_periods <- function(variable, periodos, n_units, period_length) {
  is_temp_like <- any(grepl("tavg|tmin|tmax", colnames(variable), ignore.case = TRUE))
  num_periods_calculated <- length(periodos)
  
  if (is_temp_like) {
    pnames_base <- paste0("tempe_period_", seq_len(num_periods_calculated))
  } else {
    pnames_base <- paste0("precp_period_", seq_len(num_periods_calculated))
  }
  pnames <- c(pnames_base, "min_idx", "max_idx")
  vperiods <- sapply(periodos, function(p_indices) {
    varval <- variable[, p_indices, drop = FALSE]
    if (is_temp_like) {
      return(Rfast::rowmeans(varval))
    } else {
      return(Rfast::rowsums(varval))
    }
    
  })
  if (is.vector(vperiods)) vperiods <- t(vperiods) # Handle single-row input
      
  # Find indices relative to the columns of vperiods (1 to num_periods_calculated)
  pmin <- Rfast::rowMins(vperiods, value = FALSE)
  pmax <- Rfast::rowMaxs(vperiods, value = FALSE)
  vperiods <- cbind(vperiods, pmin)
  vperiods <- cbind(vperiods, pmax)
  colnames(vperiods) <- pnames
  return(vperiods)
}


# --- Individual Bioclimatic Variable Functions (BIO 1-19) ---

#' @title BIO1: Mean Temperature of Units
#' @description Calculates mean temperature across all temporal units.
#' @param tavg Matrix of average temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with columns: "bio1", "cell".
#' @keywords internal
bio1_fun <- function(tavg, cell){
  bio1V <- Rfast::rowmeans(tavg)
  bio1V <- cbind(bio1 = bio1V, cell = cell)
  return(bio1V)
}

#' @title BIO2: Mean Diurnal Range
#' @description Calculates the mean of (tmax - tmin) across all temporal units.
#' @param tmin Matrix of minimum temperatures for each unit.
#' @param tmax Matrix of maximum temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio2", "cell".
#' @keywords internal
bio2_fun <- function(tmin, tmax, cell){
  bio2V <- tmax - tmin
  bio2V <- Rfast::rowmeans(bio2V)
  bio2V <- cbind(bio2 = bio2V, cell = cell)
  return(bio2V)
}

#' @title BIO3: Isothermality
#' @description Calculates (BIO2 / BIO7) * 100.
#' @param bio2V Vector or single-column matrix of BIO2 values.
#' @param bio7V Vector or single-column matrix of BIO7 values.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio3", "cell".
#' @keywords internal
bio3_fun <- function(bio2V, bio7V, cell){
  # Ensure inputs are treated as vectors if passed as matrices
  bio2_vec <- if(is.matrix(bio2V)) bio2V[,1] else bio2V
  bio7_vec <- if(is.matrix(bio7V)) bio7V[,1] else bio7V
  # Avoid division by zero
  bio3V <- ifelse(bio7_vec == 0, 0, (bio2_vec / bio7_vec) * 100)
  bio3V <- cbind(bio3 = bio3V, cell = cell)
  return(bio3V)
}

#' @title BIO4: Temperature Seasonality (Std Dev * 100)
#' @description Calculates the standard deviation of average temperatures across units, multiplied by 100.
#' @param tavg Matrix of average temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio4", "cell".
#' @keywords internal
bio4_fun <- function(tavg, cell){
  # std=TRUE gives standard deviation
  bio4V <- Rfast::rowVars(tavg, std = TRUE) * 100
  bio4V <- cbind(bio4 = bio4V, cell = cell)
  return(bio4V)
}

#' @title BIO5: Max Temperature of Warmest Unit
#' @description Identifies max temperature of the warmest unit, potentially using a static index.
#' @param tmax Matrix of maximum temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Tmax for that unit. If NULL, finds overall max Tmax.
#' @return Matrix with "bio5", "cell".
#' @keywords internal
bio5_fun <- function(tmax, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
    if (length(index_vector) != nrow(tmax)) stop("BIO5: Length mismatch: index_vector vs tmax rows.")
    # Handle potential NA indices or out-of-bounds indices gracefully
    is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(tmax)
    valid_rows <- which(!is_invalid_idx)
    bio5V <- rep(NA_real_, nrow(tmax)) # Initialize with NA
    if(length(valid_rows) > 0) {
        bio5V[valid_rows] <- tmax[cbind(valid_rows, index_vector[valid_rows])]
    }
    if(any(is_invalid_idx)) warning("BIO5: Some static indices were NA or out of bounds.")
  } else {
    bio5V <- Rfast::rowMaxs(tmax, value = TRUE)
  }
  bio5V <- cbind(bio5 = bio5V, cell = cell)
  return(bio5V)
}

#' @title BIO6: Min Temperature of Coldest Unit
#' @description Identifies min temperature of the coldest unit, potentially using a static index.
#' @param tmin Matrix of minimum temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Tmin for that unit. If NULL, finds overall min Tmin.
#' @return Matrix with "bio6", "cell".
#' @keywords internal
bio6_fun <- function(tmin, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(tmin)) stop("BIO6: Length mismatch: index_vector vs tmin rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(tmin)
     valid_rows <- which(!is_invalid_idx)
     bio6V <- rep(NA_real_, nrow(tmin))
     if(length(valid_rows) > 0) {
        bio6V[valid_rows] <- tmin[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("BIO6: Some static indices were NA or out of bounds.")
  } else {
    bio6V <- Rfast::rowMins(tmin, value = TRUE)
  }
  bio6V <- cbind(bio6 = bio6V, cell = cell)
  return(bio6V)
}

#' @title BIO7: Temperature Annual Range (BIO5 - BIO6)
#' @description Calculates the difference between BIO5 and BIO6.
#' @param bio5V Vector or single-column matrix of BIO5 values.
#' @param bio6V Vector or single-column matrix of BIO6 values.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio7", "cell".
#' @keywords internal
bio7_fun <- function(bio5V, bio6V, cell){
  bio5_vec <- if(is.matrix(bio5V)) bio5V[,1] else bio5V
  bio6_vec <- if(is.matrix(bio6V)) bio6V[,1] else bio6V
  bio7V <- bio5_vec - bio6_vec
  bio7V <- cbind(bio7 = bio7V, cell = cell)
  return(bio7V)
}

#' @title BIO8: Mean Temperature of Wettest Period
#' @description Calculates mean temperature of the period with the highest precipitation sum.
#' @param tperiod Matrix of temperature period sums (output from `var_periods`).
#' @param pperiod_max_idx Vector indicating the index (1-based) of the wettest period for each row.
#' @param period_length Integer. Number of units per period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio8", "cell".
#' @keywords internal
bio8_fun <- function(tperiod, pperiod_max_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2 # Exclude min/max ID columns
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO8: Some max_prec_period indices are out of bounds.")
    # Handle invalid indices - results will be NA due to matrix indexing rules
  }
  # Extract the temperature sum for the wettest period
  wettest_period_temp_sum <- tperiod[cbind(seq_len(nrow(tperiod)), pperiod_max_idx)]
  # Calculate mean
  bio8V <- wettest_period_temp_sum / period_length
  bio8V <- cbind(bio8 = bio8V, cell = cell)
  return(bio8V)
}

#' @title BIO9: Mean Temperature of Driest Period
#' @description Calculates mean temperature of the period with the lowest precipitation sum.
#' @param tperiod Matrix of temperature period sums.
#' @param pperiod_min_idx Vector indicating the index (1-based) of the driest period.
#' @param period_length Integer. Number of units per period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio9", "cell".
#' @keywords internal
bio9_fun <- function(tperiod, pperiod_min_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO9: Some min_prec_period indices are out of bounds.")
  }
  driest_period_temp_sum <- tperiod[cbind(seq_len(nrow(tperiod)), pperiod_min_idx)]
  bio9V <- driest_period_temp_sum / period_length
  bio9V <- cbind(bio9 = bio9V, cell = cell)
  return(bio9V)
}

#' @title BIO10: Mean Temperature of Warmest Period
#' @description Calculates mean temperature of the period with the highest temperature sum.
#' @param tperiod Matrix of temperature period sums.
#' @param tperiod_max_idx Vector indicating the index (1-based) of the warmest period.
#' @param period_length Integer. Number of units per period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio10", "cell".
#' @keywords internal
bio10_fun <- function(tperiod, tperiod_max_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO10: Some max_temp_period indices are out of bounds.")
  }
  warmest_period_temp_sum <- tperiod[cbind(seq_len(nrow(tperiod)), tperiod_max_idx)]
  bio10V <- warmest_period_temp_sum / period_length
  bio10V <- cbind(bio10 = bio10V, cell = cell)
  return(bio10V)
}

#' @title BIO11: Mean Temperature of Coldest Period
#' @description Calculates mean temperature of the period with the lowest temperature sum.
#' @param tperiod Matrix of temperature period sums.
#' @param tperiod_min_idx Vector indicating the index (1-based) of the coldest period.
#' @param period_length Integer. Number of units per period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio11", "cell".
#' @keywords internal
bio11_fun <- function(tperiod, tperiod_min_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO11: Some min_temp_period indices are out of bounds.")
  }
  coldest_period_temp_sum <- tperiod[cbind(seq_len(nrow(tperiod)), tperiod_min_idx)]
  bio11V <- coldest_period_temp_sum / period_length
  bio11V <- cbind(bio11 = bio11V, cell = cell)
  return(bio11V)
}

#' @title BIO12: Total Precipitation
#' @description Calculates the sum of precipitation values across all units.
#' @param precp Matrix of precipitation values for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio12", "cell".
#' @keywords internal
bio12_fun <- function(precp, cell){
  bio12V <- Rfast::rowsums(precp)
  bio12V <- cbind(bio12 = bio12V, cell = cell)
  return(bio12V)
}

#' @title BIO13: Precipitation of Wettest Unit
#' @description Identifies precipitation of the wettest unit, potentially using a static index.
#' @param precp Matrix of precipitation values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Prec for that unit. If NULL, finds overall max Prec.
#' @return Matrix with "bio13", "cell".
#' @keywords internal
bio13_fun <- function(precp, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(precp)) stop("BIO13: Length mismatch: index_vector vs precp rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(precp)
     valid_rows <- which(!is_invalid_idx)
     bio13V <- rep(NA_real_, nrow(precp))
     if(length(valid_rows) > 0) {
        bio13V[valid_rows] <- precp[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("BIO13: Some static indices were NA or out of bounds.")
  } else {
    bio13V <- Rfast::rowMaxs(precp, value = TRUE)
  }
  bio13V <- cbind(bio13 = bio13V, cell = cell)
  return(bio13V)
}

#' @title BIO14: Precipitation of Driest Unit
#' @description Identifies precipitation of the driest unit, potentially using a static index.
#' @param precp Matrix of precipitation values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Prec for that unit. If NULL, finds overall min Prec.
#' @return Matrix with "bio14", "cell".
#' @keywords internal
bio14_fun <- function(precp, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(precp)) stop("BIO14: Length mismatch: index_vector vs precp rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(precp)
     valid_rows <- which(!is_invalid_idx)
     bio14V <- rep(NA_real_, nrow(precp))
     if(length(valid_rows) > 0) {
        bio14V[valid_rows] <- precp[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("BIO14: Some static indices were NA or out of bounds.")
  } else {
    bio14V <- Rfast::rowMins(precp, value = TRUE)
  }
  bio14V <- cbind(bio14 = bio14V, cell = cell)
  return(bio14V)
}

#' @title BIO15: Precipitation Seasonality (CV)
#' @description Calculates coefficient of variation in precipitation across units.
#' @param precp Matrix containing precipitation values for each unit.
#' @param bio12V Precomputed total precipitation (BIO12 value).
#' @param n_units Integer. The total number of temporal units.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio15", "cell".
#' @keywords internal
bio15_fun <- function(precp, bio12V, n_units, cell) {
  bio12_vec <- if(is.matrix(bio12V)) bio12V[,1] else bio12V
  # Calculate mean precipitation per unit
  mean_unit_prec <- 1 + (bio12_vec / n_units) # Add 1 to total to avoid div by zero
  sd_prec <- Rfast::rowVars(precp, std = TRUE)
  # Calculate CV * 100
  bio15V <- ifelse(mean_unit_prec <= 0, 0, (sd_prec / mean_unit_prec) * 100)
  bio15V <- cbind(bio15 = bio15V, cell = cell)
  return(bio15V)
}

#' @title BIO16: Precipitation of Wettest Period
#' @description Calculates precipitation sum of the period with the highest precipitation sum.
#' @param pperiod Matrix of precipitation period sums (output from `var_periods`).
#' @param pperiod_max_idx Vector indicating the index (1-based) of the wettest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio16", "cell".
#' @keywords internal
bio16_fun <- function(pperiod, pperiod_max_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO16: Some max_prec_period indices are out of bounds.")
  }
  # Extract the sum for the wettest period
  bio16V <- pperiod[cbind(seq_len(nrow(pperiod)), pperiod_max_idx)]
  bio16V <- cbind(bio16 = bio16V, cell = cell)
  return(bio16V)
}

#' @title BIO17: Precipitation of Driest Period
#' @description Calculates precipitation sum of the period with the lowest precipitation sum.
#' @param pperiod Matrix of precipitation period sums.
#' @param pperiod_min_idx Vector indicating the index (1-based) of the driest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio17", "cell".
#' @keywords internal
bio17_fun <- function(pperiod, pperiod_min_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO17: Some min_prec_period indices are out of bounds.")
  }
  bio17V <- pperiod[cbind(seq_len(nrow(pperiod)), pperiod_min_idx)]
  bio17V <- cbind(bio17 = bio17V, cell = cell)
  return(bio17V)
}

#' @title BIO18: Precipitation of Warmest Period
#' @description Calculates precipitation sum of the period with the highest temperature sum.
#' @param pperiod Matrix of precipitation period sums.
#' @param tperiod_max_idx Vector indicating the index (1-based) of the warmest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio18", "cell".
#' @keywords internal
bio18_fun <- function(pperiod, tperiod_max_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO18: Some max_temp_period indices are out of bounds.")
  }
  bio18V <- pperiod[cbind(seq_len(nrow(pperiod)), tperiod_max_idx)]
  bio18V <- cbind(bio18 = bio18V, cell = cell)
  return(bio18V)
}

#' @title BIO19: Precipitation of Coldest Period
#' @description Calculates precipitation sum of the period with the lowest temperature sum.
#' @param pperiod Matrix of precipitation period sums.
#' @param tperiod_min_idx Vector indicating the index (1-based) of the coldest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio19", "cell".
#' @keywords internal
bio19_fun <- function(pperiod, tperiod_min_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("BIO19: Some min_temp_period indices are out of bounds.")
  }
  bio19V <- pperiod[cbind(seq_len(nrow(pperiod)), tperiod_min_idx)]
  bio19V <- cbind(bio19 = bio19V, cell = cell)
  return(bio19V)
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