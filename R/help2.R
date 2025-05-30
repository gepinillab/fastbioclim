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


# --- Individual Bioclimatic Variable Functions (bio 1-19) ---

#' @title bio01: Mean Temperature of Units
#' @description Calculates mean temperature across all temporal units.
#' @param tavg Matrix of average temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with columns: "bio01", "cell".
#' @keywords internal
bio01_fun <- function(tavg, cell){
  bio01V <- Rfast::rowmeans(tavg)
  bio01V <- cbind(bio01 = bio01V, cell = cell)
  return(bio01V)
}

#' @title bio02: Mean Diurnal Range
#' @description Calculates the mean of (tmax - tmin) across all temporal units.
#' @param tmin Matrix of minimum temperatures for each unit.
#' @param tmax Matrix of maximum temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio02", "cell".
#' @keywords internal
bio02_fun <- function(tmin, tmax, cell){
  bio02V <- tmax - tmin
  bio02V <- Rfast::rowmeans(bio02V)
  bio02V <- cbind(bio02 = bio02V, cell = cell)
  return(bio02V)
}

#' @title bio03: Isothermality
#' @description Calculates (bio02 / bio07) * 100.
#' @param bio02V Vector or single-column matrix of bio02 values.
#' @param bio07V Vector or single-column matrix of bio07 values.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio03", "cell".
#' @keywords internal
bio03_fun <- function(bio02V, bio07V, cell){
  # Ensure inputs are treated as vectors if passed as matrices
  bio02_vec <- if(is.matrix(bio02V)) bio02V[,1] else bio02V
  bio07_vec <- if(is.matrix(bio07V)) bio07V[,1] else bio07V
  # Avoid division by zero
  bio03V <- ifelse(bio07_vec == 0, 0, (bio02_vec / bio07_vec) * 100)
  bio03V <- cbind(bio03 = bio03V, cell = cell)
  return(bio03V)
}

#' @title bio04: Temperature Seasonality (Std Dev * 100)
#' @description Calculates the standard deviation of average temperatures across units, multiplied by 100.
#' @param tavg Matrix of average temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio04", "cell".
#' @keywords internal
bio04_fun <- function(tavg, cell){
  # std=TRUE gives standard deviation
  bio04V <- Rfast::rowVars(tavg, std = TRUE) * 100
  bio04V <- cbind(bio04 = bio04V, cell = cell)
  return(bio04V)
}

#' @title bio05: Max Temperature of Warmest Unit
#' @description Identifies max temperature of the warmest unit, potentially using a static index.
#' @param tmax Matrix of maximum temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Tmax for that unit. If NULL, finds overall max Tmax.
#' @return Matrix with "bio05", "cell".
#' @keywords internal
bio05_fun <- function(tmax, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
    if (length(index_vector) != nrow(tmax)) stop("bio05: Length mismatch: index_vector vs tmax rows.")
    # Handle potential NA indices or out-of-bounds indices gracefully
    is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(tmax)
    valid_rows <- which(!is_invalid_idx)
    bio05V <- rep(NA_real_, nrow(tmax)) # Initialize with NA
    if(length(valid_rows) > 0) {
        bio05V[valid_rows] <- tmax[cbind(valid_rows, index_vector[valid_rows])]
    }
    if(any(is_invalid_idx)) warning("bio05: Some static indices were NA or out of bounds.")
  } else {
    bio05V <- Rfast::rowMaxs(tmax, value = TRUE)
  }
  bio05V <- cbind(bio05 = bio05V, cell = cell)
  return(bio05V)
}

#' @title bio06: Min Temperature of Coldest Unit
#' @description Identifies min temperature of the coldest unit, potentially using a static index.
#' @param tmin Matrix of minimum temperatures for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Tmin for that unit. If NULL, finds overall min Tmin.
#' @return Matrix with "bio06", "cell".
#' @keywords internal
bio06_fun <- function(tmin, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(tmin)) stop("bio06: Length mismatch: index_vector vs tmin rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(tmin)
     valid_rows <- which(!is_invalid_idx)
     bio06V <- rep(NA_real_, nrow(tmin))
     if(length(valid_rows) > 0) {
        bio06V[valid_rows] <- tmin[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio06: Some static indices were NA or out of bounds.")
  } else {
    bio06V <- Rfast::rowMins(tmin, value = TRUE)
  }
  bio06V <- cbind(bio06 = bio06V, cell = cell)
  return(bio06V)
}

#' @title bio07: Temperature Annual Range (bio05 - bio06)
#' @description Calculates the difference between bio05 and bio06.
#' @param bio05V Vector or single-column matrix of bio05 values.
#' @param bio06V Vector or single-column matrix of bio06 values.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio07", "cell".
#' @keywords internal
bio07_fun <- function(bio05V, bio06V, cell){
  bio05_vec <- if(is.matrix(bio05V)) bio05V[,1] else bio05V
  bio06_vec <- if(is.matrix(bio06V)) bio06V[,1] else bio06V
  bio07V <- bio05_vec - bio06_vec
  bio07V <- cbind(bio07 = bio07V, cell = cell)
  return(bio07V)
}

#' @title bio08: Mean Temperature of Wettest Period
#' @description Calculates mean temperature of the period with the highest precipitation sum.
#' @param tperiod Matrix of temperature period sums (output from `var_periods`).
#' @param pperiod_max_idx Vector indicating the index (1-based) of the wettest period for each row.
#' @param period_length Integer. Number of units per period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio08", "cell".
#' @keywords internal
bio08_fun <- function(tperiod, pperiod_max_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2 # Exclude min/max ID columns
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio08: Some max_prec_period indices are out of bounds.")
    # Handle invalid indices - results will be NA due to matrix indexing rules
  }
  # Extract the temperature mean for the wettest period
  bio08V <- tperiod[cbind(seq_len(nrow(tperiod)), pperiod_max_idx)]
  bio08V <- cbind(bio08 = bio08V, cell = cell)
  return(bio08V)
}

#' @title bio09: Mean Temperature of Driest Period
#' @description Calculates mean temperature of the period with the lowest precipitation sum.
#' @param tperiod Matrix of temperature period sums.
#' @param pperiod_min_idx Vector indicating the index (1-based) of the driest period.
#' @param period_length Integer. Number of units per period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio09", "cell".
#' @keywords internal
bio09_fun <- function(tperiod, pperiod_min_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio09: Some min_prec_period indices are out of bounds.")
  }
  bio09V <- tperiod[cbind(seq_len(nrow(tperiod)), pperiod_min_idx)]
  bio09V <- cbind(bio09 = bio09V, cell = cell)
  return(bio09V)
}

#' @title bio10: Mean Temperature of Warmest Period
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
    warning("bio10: Some max_temp_period indices are out of bounds.")
  }
  bio10V <- tperiod[cbind(seq_len(nrow(tperiod)), tperiod_max_idx)]
  bio10V <- cbind(bio10 = bio10V, cell = cell)
  return(bio10V)
}

#' @title bio11: Mean Temperature of Coldest Period
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
    warning("bio11: Some min_temp_period indices are out of bounds.")
  }
  bio11V <- tperiod[cbind(seq_len(nrow(tperiod)), tperiod_min_idx)]
  bio11V <- cbind(bio11 = bio11V, cell = cell)
  return(bio11V)
}

#' @title bio12: Total Precipitation
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

#' @title bio13: Precipitation of Wettest Unit
#' @description Identifies precipitation of the wettest unit, potentially using a static index.
#' @param precp Matrix of precipitation values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Prec for that unit. If NULL, finds overall max Prec.
#' @return Matrix with "bio13", "cell".
#' @keywords internal
bio13_fun <- function(precp, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(precp)) stop("bio13: Length mismatch: index_vector vs precp rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(precp)
     valid_rows <- which(!is_invalid_idx)
     bio13V <- rep(NA_real_, nrow(precp))
     if(length(valid_rows) > 0) {
        bio13V[valid_rows] <- precp[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio13: Some static indices were NA or out of bounds.")
  } else {
    bio13V <- Rfast::rowMaxs(precp, value = TRUE)
  }
  bio13V <- cbind(bio13 = bio13V, cell = cell)
  return(bio13V)
}

#' @title bio14: Precipitation of Driest Unit
#' @description Identifies precipitation of the driest unit, potentially using a static index.
#' @param precp Matrix of precipitation values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Prec for that unit. If NULL, finds overall min Prec.
#' @return Matrix with "bio14", "cell".
#' @keywords internal
bio14_fun <- function(precp, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(precp)) stop("bio14: Length mismatch: index_vector vs precp rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(precp)
     valid_rows <- which(!is_invalid_idx)
     bio14V <- rep(NA_real_, nrow(precp))
     if(length(valid_rows) > 0) {
        bio14V[valid_rows] <- precp[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio14: Some static indices were NA or out of bounds.")
  } else {
    bio14V <- Rfast::rowMins(precp, value = TRUE)
  }
  bio14V <- cbind(bio14 = bio14V, cell = cell)
  return(bio14V)
}

#' @title bio15: Precipitation Seasonality (CV)
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

#' @title bio16: Precipitation of Wettest Period
#' @description Calculates precipitation sum of the period with the highest precipitation sum.
#' @param pperiod Matrix of precipitation period sums (output from `var_periods`).
#' @param pperiod_max_idx Vector indicating the index (1-based) of the wettest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio16", "cell".
#' @keywords internal
bio16_fun <- function(pperiod, pperiod_max_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio16: Some max_prec_period indices are out of bounds.")
  }
  # Extract the sum for the wettest period
  bio16V <- pperiod[cbind(seq_len(nrow(pperiod)), pperiod_max_idx)]
  bio16V <- cbind(bio16 = bio16V, cell = cell)
  return(bio16V)
}

#' @title bio17: Precipitation of Driest Period
#' @description Calculates precipitation sum of the period with the lowest precipitation sum.
#' @param pperiod Matrix of precipitation period sums.
#' @param pperiod_min_idx Vector indicating the index (1-based) of the driest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio17", "cell".
#' @keywords internal
bio17_fun <- function(pperiod, pperiod_min_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio17: Some min_prec_period indices are out of bounds.")
  }
  bio17V <- pperiod[cbind(seq_len(nrow(pperiod)), pperiod_min_idx)]
  bio17V <- cbind(bio17 = bio17V, cell = cell)
  return(bio17V)
}

#' @title bio18: Precipitation of Warmest Period
#' @description Calculates precipitation sum of the period with the highest temperature sum.
#' @param pperiod Matrix of precipitation period sums.
#' @param tperiod_max_idx Vector indicating the index (1-based) of the warmest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio18", "cell".
#' @keywords internal
bio18_fun <- function(pperiod, tperiod_max_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio18: Some max_temp_period indices are out of bounds.")
  }
  bio18V <- pperiod[cbind(seq_len(nrow(pperiod)), tperiod_max_idx)]
  bio18V <- cbind(bio18 = bio18V, cell = cell)
  return(bio18V)
}

#' @title bio19: Precipitation of Coldest Period
#' @description Calculates precipitation sum of the period with the lowest temperature sum.
#' @param pperiod Matrix of precipitation period sums.
#' @param tperiod_min_idx Vector indicating the index (1-based) of the coldest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio19", "cell".
#' @keywords internal
bio19_fun <- function(pperiod, tperiod_min_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio19: Some min_temp_period indices are out of bounds.")
  }
  bio19V <- pperiod[cbind(seq_len(nrow(pperiod)), tperiod_min_idx)]
  bio19V <- cbind(bio19 = bio19V, cell = cell)
  return(bio19V)
}

#' @title bio20: Mean Solar Radiation of Units
#' @description Calculates mean solar radiation across all temporal units.
#' @param srad Matrix of average solar radiation for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with columns: "bio20", "cell".
#' @keywords internal
bio20_fun <- function(srad, cell) {
  bio20V <- Rfast::rowmeans(srad)
  bio20V <- cbind(bio20 = bio20V, cell = cell)
  return(bio20V)
}

#' @title bio21: Highest Solar Radiation Unit
#' @description Identifies highest solar radiation unit, potentially using a static index.
#' @param srad Matrix of solar radiation values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Srad for that unit. If NULL, finds overall max Srad.
#' @return Matrix with "bio21", "cell".
#' @keywords internal
bio21_fun <- function(srad, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(srad)) stop("bio21: Length mismatch: index_vector vs srad rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(srad)
     valid_rows <- which(!is_invalid_idx)
     bio21V <- rep(NA_real_, nrow(srad))
     if(length(valid_rows) > 0) {
        bio21V[valid_rows] <- srad[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio21: Some static indices were NA or out of bounds.")
  } else {
    bio21V <- Rfast::rowMaxs(srad, value = TRUE)
  }
  bio21V <- cbind(bio21 = bio21V, cell = cell)
  return(bio21V)
}

#' @title bio22: Lowest Solar Radiation Unit
#' @description Identifies lowest solar radiation unit, potentially using a static index.
#' @param srad Matrix of solar radiation values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts Srad for that unit. If NULL, finds overall max Srad.
#' @return Matrix with "bio22", "cell".
#' @keywords internal
bio22_fun <- function(srad, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(srad)) stop("bio22: Length mismatch: index_vector vs srad rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(srad)
     valid_rows <- which(!is_invalid_idx)
     bio22V <- rep(NA_real_, nrow(srad))
     if(length(valid_rows) > 0) {
        bio22V[valid_rows] <- srad[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio22: Some static indices were NA or out of bounds.")
  } else {
    bio22V <- Rfast::rowMins(srad, value = TRUE)
  }
  bio22V <- cbind(bio22 = bio22V, cell = cell)
  return(bio22V)
}

#' @title bio23: Solar Radiation Seasonality (CV)
#' @description Calculates coefficient of variation in solar radiation across units.
#' @param srad Matrix containing solar radiation values for each unit.
#' @param n_units Integer. The total number of temporal units.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio23", "cell".
#' @keywords internal
bio23_fun <- function(srad, n_units, cell) {
  # Calculate mean solar radiation
  mean_unit_srad <- 1 + Rfast::rowmeans(srad) # Add 1 to total to avoid div by zero
  sd_srad <- Rfast::rowVars(srad, std = TRUE)
  # Calculate CV * 100
  bio23V <- ifelse(mean_unit_srad <= 0, 0, (sd_srad / mean_unit_srad) * 100)
  bio23V <- cbind(bio23 = bio23V, cell = cell)
  return(bio23V)
}

#' @title bio24: Solar Radiation of Wettest Period
#' @description Calculates solar radiation mean of the period with the highest precipitation sum.
#' @param speriod Matrix of solar radiation period means (output from `var_periods`).
#' @param pperiod_max_idx Vector indicating the index (1-based) of the wettest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio24", "cell".
#' @keywords internal
bio24_fun <- function(speriod, pperiod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio24: Some max_prec_period indices are out of bounds.")
  }
  # Extract the mean for the wettest period
  bio24V <- speriod[cbind(seq_len(nrow(speriod)), pperiod_max_idx)]
  bio24V <- cbind(bio24 = bio24V, cell = cell)
  return(bio24V)
}

#' @title bio25: Solar Radiation of Driest Period
#' @description Calculates solar radiation mean of the period with the highest precipitation sum.
#' @param speriod Matrix of solar radiation period means (output from `var_periods`).
#' @param pperiod_min_idx Vector indicating the index (1-based) of the driest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio25", "cell".
#' @keywords internal
bio25_fun <- function(speriod, pperiod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio25: Some min_prec_period indices are out of bounds.")
  }
  # Extract the mean for the driest period
  bio25V <- speriod[cbind(seq_len(nrow(speriod)), pperiod_min_idx)]
  bio25V <- cbind(bio25 = bio25V, cell = cell)
  return(bio25V)
}

#' @title bio26: Solar Radiation of Warmest Period
#' @description Calculates solar radiation mean of the period with the highest temperature mean.
#' @param speriod Matrix of solar radiation period means (output from `var_periods`).
#' @param tperiod_max_idx Vector indicating the index (1-based) of the warmest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio26", "cell".
#' @keywords internal
bio26_fun <- function(speriod, tperiod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio26: Some max_temp_period indices are out of bounds.")
  }
  # Extract the mean for the warmest period
  bio26V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_max_idx)]
  bio26V <- cbind(bio26 = bio26V, cell = cell)
  return(bio26V)
}

#' @title bio27: Solar Radiation of Coldest Period
#' @description Calculates solar radiation mean of the period with the lowest temperature mean.
#' @param speriod Matrix of solar radiation period means (output from `var_periods`).
#' @param tperiod_min_idx Vector indicating the index (1-based) of the coldest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio27", "cell".
#' @keywords internal
bio27_fun <- function(speriod, tperiod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio27: Some min_temp_period indices are out of bounds.")
  }
  # Extract the mean for the coldest period
  bio27V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_min_idx)]
  bio27V <- cbind(bio27 = bio27V, cell = cell)
  return(bio27V)
}

#' @title bio28: Mean Moisture of Units
#' @description Calculates mean moisture across all temporal units.
#' @param mois Matrix of average moisture for each unit.
#' @param cell Vector of original cell IDs.
#' @return Matrix with columns: "bio28", "cell".
#' @keywords internal
bio28_fun <- function(mois, cell) {
  bio28V <- Rfast::rowmeans(mois)
  bio28V <- cbind(bio28 = bio28V, cell = cell)
  return(bio28V)
}

#' @title bio29: Highest Moisture Unit
#' @description Identifies highest moisture unit, potentially using a static index.
#' @param mois Matrix of moisture values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts mois for that unit. If NULL, finds overall max mois.
#' @return Matrix with "bio29", "cell".
#' @keywords internal
bio29_fun <- function(mois, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(mois)) stop("bio29: Length mismatch: index_vector vs mois rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(mois)
     valid_rows <- which(!is_invalid_idx)
     bio29V <- rep(NA_real_, nrow(mois))
     if(length(valid_rows) > 0) {
        bio29V[valid_rows] <- mois[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio29: Some static indices were NA or out of bounds.")
  } else {
    bio29V <- Rfast::rowMaxs(mois, value = TRUE)
  }
  bio29V <- cbind(bio29 = bio29V, cell = cell)
  return(bio29V)
}

#' @title bio30: Lowest Moisture Unit
#' @description Identifies lowest moisture unit, potentially using a static index.
#' @param mois Matrix of moisture values for each unit.
#' @param cell Vector of original cell IDs.
#' @param index_vector Optional vector of unit indices (1-based). If provided, extracts mois for that unit. If NULL, finds overall max mois.
#' @return Matrix with "bio30", "cell".
#' @keywords internal
bio30_fun <- function(mois, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(mois)) stop("bio30: Length mismatch: index_vector vs mois rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(mois)
     valid_rows <- which(!is_invalid_idx)
     bio30V <- rep(NA_real_, nrow(mois))
     if(length(valid_rows) > 0) {
        bio30V[valid_rows] <- mois[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio30: Some static indices were NA or out of bounds.")
  } else {
    bio30V <- Rfast::rowMins(mois, value = TRUE)
  }
  bio30V <- cbind(bio30 = bio30V, cell = cell)
  return(bio30V)
}

#' @title bio31: Moisture Seasonality (Standard Deviation)
#' @description Calculates coefficient of variation in moisture across units.
#' @param mois Matrix containing moisture values for each unit.
#' @param n_units Integer. The total number of temporal units.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio31", "cell".
#' @keywords internal
bio31_fun <- function(mois, n_units, cell) {
  bio31V <- Rfast::rowVars(mois, std = TRUE) * 100
  bio31V <- cbind(bio31 = bio31V, cell = cell)
  return(bio31V)
}

#' @title bio32: Moisture of the Most Moist Period
#' @description Calculates moisture mean of the most moist period.
#' @param speriod Matrix of moisture period means (output from `var_periods`).
#' @param speriod_max_idx Vector indicating the index (1-based) of the most moist period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio32", "cell".
#' @keywords internal
bio32_fun <- function(speriod, speriod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(speriod_max_idx < 1, na.rm=TRUE) || any(speriod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio32: Some max_prec_period indices are out of bounds.")
  }
  # Extract the mean for the most moist period
  bio32V <- speriod[cbind(seq_len(nrow(speriod)), speriod_max_idx)]
  bio32V <- cbind(bio32 = bio32V, cell = cell)
  return(bio32V)
}

#' @title bio33: Moisture of the Least Moist Period
#' @description Calculates moisture mean of the least moist period.
#' @param speriod Matrix of moisture period means (output from `var_periods`).
#' @param speriod_min_idx Vector indicating the index (1-based) of the least moist period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio33", "cell".
#' @keywords internal
bio33_fun <- function(speriod, speriod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(speriod_min_idx < 1, na.rm=TRUE) || any(speriod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio33: Some min_prec_period indices are out of bounds.")
  }
  # Extract the mean for the least moist period
  bio33V <- speriod[cbind(seq_len(nrow(speriod)), speriod_min_idx)]
  bio33V <- cbind(bio33 = bio33V, cell = cell)
  return(bio33V)
}

#' @title bio34: Moisture of Warmest Period
#' @description Calculates moisture mean of the period with the highest temperature mean.
#' @param speriod Matrix of moisture period means (output from `var_periods`).
#' @param tperiod_max_idx Vector indicating the index (1-based) of the warmest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio34", "cell".
#' @keywords internal
bio34_fun <- function(speriod, tperiod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio34: Some max_temp_period indices are out of bounds.")
  }
  # Extract the mean for the warmest period
  bio34V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_max_idx)]
  bio34V <- cbind(bio34 = bio34V, cell = cell)
  return(bio34V)
}

#' @title bio35: Moisture of Coldest Period
#' @description Calculates moisture mean of the period with the lowest temperature mean.
#' @param speriod Matrix of moisture period means (output from `var_periods`).
#' @param tperiod_min_idx Vector indicating the index (1-based) of the coldest period.
#' @param cell Vector of original cell IDs.
#' @return Matrix with "bio35", "cell".
#' @keywords internal
bio35_fun <- function(speriod, tperiod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio35: Some min_temp_period indices are out of bounds.")
  }
  # Extract the mean for the coldest period
  bio35V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_min_idx)]
  bio35V <- cbind(bio35 = bio35V, cell = cell)
  return(bio35V)
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