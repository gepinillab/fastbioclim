#' @title bio01_fast: Mean Temperature
#' @description Calculates mean temperature across all temporal units (usually 12 months).
#' @param tavg A numeric **matrix** where **rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., 12 months).
#' @param cell A numeric or character **vector** of original cell IDs. 
#'   Its length must be exactly equal to the number of rows in `tavg`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio01" (the calculated mean) and "cell" (the IDs).
#' @keywords internal
bio01_fast <- function(tavg, cell){
  bio01V <- Rfast::rowmeans(tavg)
  bio01V <- cbind(bio01 = bio01V, cell = cell)
  return(bio01V)
}

#' @title bio02_fast: Mean Diurnal Range
#' @description Calculates the mean of the diurnal temperature ranges (tmax - tmin) across all temporal units.
#' @param tmin A numeric **matrix** of minimum temperatures. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units. Must have the exact same dimensions as `tmax`.
#' @param tmax A numeric **matrix** of maximum temperatures. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units. Must have the exact same dimensions as `tmin`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tmin` and `tmax`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio02" (mean diurnal range) and "cell".
#' @keywords internal
bio02_fast <- function(tmin, tmax, cell){
  bio02V <- tmax - tmin
  bio02V <- Rfast::rowmeans(bio02V)
  bio02V <- cbind(bio02 = bio02V, cell = cell)
  return(bio02V)
}

#' @title bio03_fast: Isothermality
#' @description Calculates Isothermality representing the ratio of mean diurnal range to temperature range: (bio02 / bio07) * 100.
#' @param bio02V A numeric **vector** or **single-column matrix** of Bio02 values (Mean Diurnal Range). 
#'   Length (or number of rows) must match `bio07V` and `cell`.
#' @param bio07V A numeric **vector** or **single-column matrix** of Bio07 values (Temperature Range). 
#'   Length (or number of rows) must match `bio02V` and `cell`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the length/rows of `bio02V` and `bio07V`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio03" and "cell".
#' @keywords internal
bio03_fast <- function(bio02V, bio07V, cell){
  # Ensure inputs are treated as vectors if passed as matrices
  bio02_vec <- if(is.matrix(bio02V)) bio02V[,1] else bio02V
  bio07_vec <- if(is.matrix(bio07V)) bio07V[,1] else bio07V
  # Avoid division by zero
  bio03V <- ifelse(bio07_vec == 0, 0, (bio02_vec / bio07_vec) * 100)
  bio03V <- cbind(bio03 = bio03V, cell = cell)
  return(bio03V)
}

#' @title bio04_fast: Temperature Seasonality (Std Dev * 100)
#' @description Calculates Temperature Seasonality, defined as the standard deviation of average temperatures across all temporal units (e.g., 12 months), multiplied by 100.
#' @param tavg A numeric **matrix** of average temperatures. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tavg`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio04" (seasonality) and "cell".
#' @keywords internal
bio04_fast <- function(tavg, cell){
  # std=TRUE gives standard deviation
  bio04V <- Rfast::rowVars(tavg, std = TRUE) * 100
  bio04V <- cbind(bio04 = bio04V, cell = cell)
  return(bio04V)
}

#' @title bio05_fast: Max Temperature of Warmest Unit
#' @description Identifies the maximum temperature of the warmest temporal unit. 
#'   If `index_vector` is `NULL`, it calculates the row-wise maximum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @details This function calculates bio05 following the standard definition used by WorldClim 
#'   (Hijmans et al., 2005) and ANUCLIM 6.1 (Xu & Hutchinson, 2013), which is the single highest value from all 
#'   maximum temperature layers (a "max of maxes"). It does not use mean 
#'   temperature to first identify the warmest month.
#' @param tmax A numeric **matrix** of maximum temperatures. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tmax`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `tmax`. 
#'   Values must be between 1 and `ncol(tmax)`. 
#'   This is typically used to extract the Tmax of the specific month identified as the warmest by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio05" (the maximum temperature) and "cell".
#' @references
#' Hijmans, R.J., Cameron, S.E., Parra, J.L., Jones, P.G. and Jarvis, A. (2005). 
#' Very high resolution interpolated climate surfaces for global land areas. 
#' International Journal of Climatology, 25(15), 1965-1978.
#'
#' O'Donnell, M. S., & Ignizio, D. A. (2012). Bioclimatic predictors for 
#' supporting ecological applications in the conterminous United States. 
#' U.S. Geological Survey Data Series 691.
#'
#' Xu, T., & Hutchinson, M. F. (2013). New developments and applications in the 
#' ANUCLIM spatial climatic and bioclimatic modelling package. Environmental 
#' Modelling & Software, 40, 267-279.
#' @keywords internal
bio05_fast <- function(tmax, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
    if (length(index_vector) != nrow(tmax)) stop("bio05_fast: Length mismatch: index_vector vs tmax rows.")
    # Handle potential NA indices or out-of-bounds indices gracefully
    is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(tmax)
    valid_rows <- which(!is_invalid_idx)
    bio05V <- rep(NA_real_, nrow(tmax)) # Initialize with NA
    if(length(valid_rows) > 0) {
        bio05V[valid_rows] <- tmax[cbind(valid_rows, index_vector[valid_rows])]
    }
    if(any(is_invalid_idx)) warning("bio05_fast: Some static indices were NA or out of bounds.")
  } else {
    bio05V <- Rfast::rowMaxs(tmax, value = TRUE)
  }
  bio05V <- cbind(bio05 = bio05V, cell = cell)
  return(bio05V)
}

#' @title bio06_fast: Min Temperature of Coldest Unit
#' @description Identifies the minimum temperature of the coldest temporal unit. 
#'   If `index_vector` is `NULL`, it calculates the row-wise minimum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @details This function calculates bio06 following the standard definition used by WorldClim 
#'   (Hijmans et al., 2005) and ANUCLIM 6.1 (Xu & Hutchinson, 2013), which is the single lowest value from all 
#'   minimum temperature layers (a "min of mins"). It does not use mean 
#'   temperature to first identify the coldest month.
#' @param tmin A numeric **matrix** of minimum temperatures. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tmin`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `tmin`. 
#'   Values must be between 1 and `ncol(tmin)`. 
#'   This is typically used to extract the Tmin of the specific month identified as the coldest by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio06" (the minimum temperature) and "cell".
#' @references
#' Hijmans, R.J., Cameron, S.E., Parra, J.L., Jones, P.G. and Jarvis, A. (2005). 
#' Very high resolution interpolated climate surfaces for global land areas. 
#' International Journal of Climatology, 25(15), 1965-1978.
#'
#' O'Donnell, M. S., & Ignizio, D. A. (2012). Bioclimatic predictors for 
#' supporting ecological applications in the conterminous United States. 
#' U.S. Geological Survey Data Series 691.
#'
#' Xu, T., & Hutchinson, M. F. (2013). New developments and applications in the 
#' ANUCLIM spatial climatic and bioclimatic modelling package. Environmental 
#' Modelling & Software, 40, 267-279.
#' @keywords internal
bio06_fast <- function(tmin, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(tmin)) stop("bio06_fast: Length mismatch: index_vector vs tmin rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(tmin)
     valid_rows <- which(!is_invalid_idx)
     bio06V <- rep(NA_real_, nrow(tmin))
     if(length(valid_rows) > 0) {
        bio06V[valid_rows] <- tmin[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio06_fast: Some static indices were NA or out of bounds.")
  } else {
    bio06V <- Rfast::rowMins(tmin, value = TRUE)
  }
  bio06V <- cbind(bio06 = bio06V, cell = cell)
  return(bio06V)
}

#' @title bio07_fast: Temperature Range (bio05 - bio06)
#' @description Calculates the Temperature Range, defined as the difference between the Maximum Temperature of the Warmest Unit (bio05) and the Minimum Temperature of the Coldest Unit (bio06).
#' @param bio05V A numeric **vector** or **single-column matrix** of Bio05 values. 
#'   Length (or number of rows) must match `bio06V` and `cell`.
#' @param bio06V A numeric **vector** or **single-column matrix** of Bio06 values. 
#'   Length (or number of rows) must match `bio05V` and `cell`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the length/rows of `bio05V` and `bio06V`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio07" (annual range) and "cell".
#' @keywords internal
bio07_fast <- function(bio05V, bio06V, cell){
  bio05_vec <- if(is.matrix(bio05V)) bio05V[,1] else bio05V
  bio06_vec <- if(is.matrix(bio06V)) bio06V[,1] else bio06V
  bio07V <- bio05_vec - bio06_vec
  bio07V <- cbind(bio07 = bio07V, cell = cell)
  return(bio07V)
}

#' @title bio08_fast: Mean Temperature of Wettest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the wettest (highest precipitation).
#' @param tperiod A numeric **matrix** of temperature values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods (typically 12). 
#' @param pperiod_max_idx An integer **vector** indicating the column index (1-based) of the wettest period for each row. 
#'   Its length must be exactly equal to the number of rows in `tperiod`.
#' @param period_length A single **integer** representing the number of units per period (e.g., 3 for months). 
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio08" (mean temperature of wettest period) and "cell".
#' @keywords internal
bio08_fast <- function(tperiod, pperiod_max_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2 # Exclude min/max ID columns
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio08_fast: Some max_prec_period indices are out of bounds.")
    # Handle invalid indices - results will be NA due to matrix indexing rules
  }
  # Extract the temperature mean for the wettest period
  bio08V <- tperiod[cbind(seq_len(nrow(tperiod)), pperiod_max_idx)]
  bio08V <- cbind(bio08 = bio08V, cell = cell)
  return(bio08V)
}

#' @title bio09_fast: Mean Temperature of Driest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the driest (lowest precipitation).
#' @param tperiod A numeric **matrix** of temperature values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods.
#' @param pperiod_min_idx An integer **vector** indicating the column index (1-based) of the driest period for each row. 
#'   Its length must be exactly equal to the number of rows in `tperiod`.
#' @param period_length A single **integer** representing the number of units per period. 
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio09" (mean temperature of driest period) and "cell".
#' @keywords internal
bio09_fast <- function(tperiod, pperiod_min_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio09_fast: Some min_prec_period indices are out of bounds.")
  }
  bio09V <- tperiod[cbind(seq_len(nrow(tperiod)), pperiod_min_idx)]
  bio09V <- cbind(bio09 = bio09V, cell = cell)
  return(bio09V)
}

#' @title bio10_fast: Mean Temperature of Warmest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the warmest (highest temperature).
#' @param tperiod A numeric **matrix** of temperature values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods.
#' @param tperiod_max_idx An integer **vector** indicating the column index (1-based) of the warmest period for each row. 
#'   Its length must be exactly equal to the number of rows in `tperiod`.
#' @param period_length A single **integer** representing the number of units per period. 
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio10" (mean temperature of warmest period) and "cell".
#' @keywords internal
bio10_fast <- function(tperiod, tperiod_max_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio10_fast: Some max_temp_period indices are out of bounds.")
  }
  bio10V <- tperiod[cbind(seq_len(nrow(tperiod)), tperiod_max_idx)]
  bio10V <- cbind(bio10 = bio10V, cell = cell)
  return(bio10V)
}

#' @title bio11_fast: Mean Temperature of Coldest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the coldest (lowest temperature).
#' @param tperiod A numeric **matrix** of temperature values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods.
#' @param tperiod_min_idx An integer **vector** indicating the column index (1-based) of the coldest period for each row. 
#'   Its length must be exactly equal to the number of rows in `tperiod`.
#' @param period_length A single **integer** representing the number of units per period. 
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `tperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio11" (mean temperature of coldest period) and "cell".
#' @keywords internal
bio11_fast <- function(tperiod, tperiod_min_idx, period_length, cell){
  num_period_cols <- ncol(tperiod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio11_fast: Some min_temp_period indices are out of bounds.")
  }
  bio11V <- tperiod[cbind(seq_len(nrow(tperiod)), tperiod_min_idx)]
  bio11V <- cbind(bio11 = bio11V, cell = cell)
  return(bio11V)
}

#' @title bio12_fast: Total Precipitation
#' @description Calculates the total precipitation (sum) across all temporal units (e.g., 12 months).
#' @param prcp A numeric **matrix** of precipitation values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `prcp`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio12" (total precipitation) and "cell".
#' 
#' @keywords internal
bio12_fast <- function(prcp, cell){
  bio12V <- Rfast::rowsums(prcp)
  bio12V <- cbind(bio12 = bio12V, cell = cell)
  return(bio12V)
}

#' @title bio13_fast: Precipitation of Wettest Unit
#' @description Identifies the precipitation of the wettest temporal unit (e.g., month). 
#'   If `index_vector` is `NULL`, it calculates the row-wise maximum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @param prcp A numeric **matrix** of precipitation values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `prcp`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `prcp`. 
#'   Values must be between 1 and `ncol(prcp)`. 
#'   This is typically used to extract the precipitation of the specific month identified as the wettest by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio13" (precipitation of wettest unit) and "cell".
#' @keywords internal
bio13_fast <- function(prcp, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(prcp)) stop("bio13_fast: Length mismatch: index_vector vs prcp rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(prcp)
     valid_rows <- which(!is_invalid_idx)
     bio13V <- rep(NA_real_, nrow(prcp))
     if(length(valid_rows) > 0) {
        bio13V[valid_rows] <- prcp[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio13_fast: Some static indices were NA or out of bounds.")
  } else {
    bio13V <- Rfast::rowMaxs(prcp, value = TRUE)
  }
  bio13V <- cbind(bio13 = bio13V, cell = cell)
  return(bio13V)
}

#' @title bio14_fast: Precipitation of Driest Unit
#' @description Identifies the precipitation of the driest temporal unit (e.g., month). 
#'   If `index_vector` is `NULL`, it calculates the row-wise minimum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @param prcp A numeric **matrix** of precipitation values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `prcp`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `prcp`. 
#'   Values must be between 1 and `ncol(prcp)`. 
#'   This is typically used to extract the precipitation of the specific month identified as the driest by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio14" (precipitation of driest unit) and "cell".
#' @keywords internal
bio14_fast <- function(prcp, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(prcp)) stop("bio14_fast: Length mismatch: index_vector vs prcp rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(prcp)
     valid_rows <- which(!is_invalid_idx)
     bio14V <- rep(NA_real_, nrow(prcp))
     if(length(valid_rows) > 0) {
        bio14V[valid_rows] <- prcp[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio14_fast: Some static indices were NA or out of bounds.")
  } else {
    bio14V <- Rfast::rowMins(prcp, value = TRUE)
  }
  bio14V <- cbind(bio14 = bio14V, cell = cell)
  return(bio14V)
}

#' @title bio15_fast: Precipitation Seasonality (CV)
#' @description Calculates the Coefficient of Variation (CV) of precipitation. 
#'   The formula used is: `(StandardDeviation / (Mean + 1)) * 100`. 
#'   (The "+1" is added to the mean to avoid division by zero in completely arid areas).
#' @param prcp A numeric **matrix** of precipitation values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., 12 months).
#' @param bio12V A numeric **vector** or **single-column matrix** of Total Annual Precipitation (Bio12). 
#'   Its length (or number of rows) must be exactly equal to the number of rows in `prcp`.
#' @param n_units A single **integer** representing the number of temporal units (e.g., 12). 
#'   This is used to calculate the mean precipitation from the total `bio12V`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `prcp`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio15" (Precipitation Seasonality) and "cell".
#' @keywords internal
bio15_fast <- function(prcp, bio12V, n_units, cell) {
  bio12_vec <- if(is.matrix(bio12V)) bio12V[,1] else bio12V
  # Calculate mean precipitation per unit
  mean_unit_prec <- 1 + (bio12_vec / n_units) # Add 1 to total to avoid div by zero
  sd_prec <- Rfast::rowVars(prcp, std = TRUE)
  # Calculate CV * 100
  bio15V <- ifelse(mean_unit_prec <= 0, 0, (sd_prec / mean_unit_prec) * 100)
  bio15V <- cbind(bio15 = bio15V, cell = cell)
  return(bio15V)
}

#' @title bio16_fast: Precipitation of Wettest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the wettest (highest precipitation).
#' @param pperiod A numeric **matrix** of precipitation sums for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param pperiod_max_idx An integer **vector** indicating the column index (1-based) of the wettest period for each row. 
#'   Its length must be exactly equal to the number of rows in `pperiod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `pperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio16" (precipitation of wettest period) and "cell".
#' @keywords internal
bio16_fast <- function(pperiod, pperiod_max_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio16_fast: Some max_prec_period indices are out of bounds.")
  }
  # Extract the sum for the wettest period
  bio16V <- pperiod[cbind(seq_len(nrow(pperiod)), pperiod_max_idx)]
  bio16V <- cbind(bio16 = bio16V, cell = cell)
  return(bio16V)
}

#' @title bio17_fast: Precipitation of Driest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the driest (lowest precipitation).
#' @param pperiod A numeric **matrix** of precipitation sums for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param pperiod_min_idx An integer **vector** indicating the column index (1-based) of the driest period for each row. 
#'   Its length must be exactly equal to the number of rows in `pperiod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `pperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio17" (precipitation of driest period) and "cell".
#' @keywords internal
bio17_fast <- function(pperiod, pperiod_min_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio17_fast: Some min_prec_period indices are out of bounds.")
  }
  bio17V <- pperiod[cbind(seq_len(nrow(pperiod)), pperiod_min_idx)]
  bio17V <- cbind(bio17 = bio17V, cell = cell)
  return(bio17V)
}

#' @title bio18_fast: Precipitation of Warmest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the warmest (highest temperature).
#' @param pperiod A numeric **matrix** of precipitation sums for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param tperiod_max_idx An integer **vector** indicating the column index (1-based) of the warmest period for each row. 
#'   Its length must be exactly equal to the number of rows in `pperiod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `pperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio18" (precipitation of warmest period) and "cell".
#' @keywords internal
bio18_fast <- function(pperiod, tperiod_max_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio18_fast: Some max_temp_period indices are out of bounds.")
  }
  bio18V <- pperiod[cbind(seq_len(nrow(pperiod)), tperiod_max_idx)]
  bio18V <- cbind(bio18 = bio18V, cell = cell)
  return(bio18V)
}

#' @title bio19_fast: Precipitation of Coldest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the coldest (lowest temperature).
#' @param pperiod A numeric **matrix** of precipitation sums for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param tperiod_min_idx An integer **vector** indicating the column index (1-based) of the coldest period for each row. 
#'   Its length must be exactly equal to the number of rows in `pperiod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `pperiod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio19" (precipitation of coldest period) and "cell".
#' @keywords internal
bio19_fast <- function(pperiod, tperiod_min_idx, cell){
  num_period_cols <- ncol(pperiod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio19_fast: Some min_temp_period indices are out of bounds.")
  }
  bio19V <- pperiod[cbind(seq_len(nrow(pperiod)), tperiod_min_idx)]
  bio19V <- cbind(bio19 = bio19V, cell = cell)
  return(bio19V)
}

#' @title bio20_fast: Mean Radiation
#' @description Calculates mean solar radiation across all temporal units (usually 12 months).
#' @param srad A numeric **matrix** where **rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., 12 months).
#' @param cell A numeric or character **vector** of original cell IDs. 
#'   Its length must be exactly equal to the number of rows in `srad`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio20" (the calculated mean) and "cell" (the IDs).
#' @keywords internal
bio20_fast <- function(srad, cell) {
  bio20V <- Rfast::rowmeans(srad)
  bio20V <- cbind(bio20 = bio20V, cell = cell)
  return(bio20V)
}

#' #' @title bio21_fast: Highest Radiation Unit
#' @description Identifies the highest solar radiation of the temporal unit with the highest value. 
#'   If `index_vector` is `NULL`, it calculates the row-wise maximum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @param srad A numeric **matrix** of solar radiation values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `srad`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `srad`. 
#'   Values must be between 1 and `ncol(srad)`. 
#'   This is typically used to extract the solar radiation of a specific unit identified by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio21" (the maximum solar radiation) and "cell".
#' @keywords internal
bio21_fast <- function(srad, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(srad)) stop("bio21_fast: Length mismatch: index_vector vs srad rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(srad)
     valid_rows <- which(!is_invalid_idx)
     bio21V <- rep(NA_real_, nrow(srad))
     if(length(valid_rows) > 0) {
        bio21V[valid_rows] <- srad[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio21_fast: Some static indices were NA or out of bounds.")
  } else {
    bio21V <- Rfast::rowMaxs(srad, value = TRUE)
  }
  bio21V <- cbind(bio21 = bio21V, cell = cell)
  return(bio21V)
}

#' @title bio22_fast: Lowest Radiation Unit
#' @description Identifies the lowest solar radiation of the temporal unit with the lowest value. 
#'   If `index_vector` is `NULL`, it calculates the row-wise minimum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @param srad A numeric **matrix** of solar radiation values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `srad`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `srad`. 
#'   Values must be between 1 and `ncol(srad)`. 
#'   This is typically used to extract the solar radiation of a specific unit identified by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio22" (the minimum solar radiation) and "cell".
#' @keywords internal
bio22_fast <- function(srad, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(srad)) stop("bio22_fast: Length mismatch: index_vector vs srad rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(srad)
     valid_rows <- which(!is_invalid_idx)
     bio22V <- rep(NA_real_, nrow(srad))
     if(length(valid_rows) > 0) {
        bio22V[valid_rows] <- srad[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio22_fast: Some static indices were NA or out of bounds.")
  } else {
    bio22V <- Rfast::rowMins(srad, value = TRUE)
  }
  bio22V <- cbind(bio22 = bio22V, cell = cell)
  return(bio22V)
}

#' @title bio23_fast: Radiation Seasonality (CV)
#' @description Calculates the Coefficient of Variation (CV) of solar radiation. 
#'   The formula used is: `(StandardDeviation / (Mean + 1)) * 100`. 
#'   (The "+1" is added to the mean to avoid division by zero).
#' @param srad A numeric **matrix** of solar radiation values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., 12 months).
#' @param n_units A single **integer** representing the number of temporal units (e.g., 12).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `srad`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio23" (Solar Radiation Seasonality) and "cell".
#' @keywords internal
bio23_fast <- function(srad, n_units, cell) {
  # Calculate mean solar radiation
  mean_unit_srad <- 1 + Rfast::rowmeans(srad) # Add 1 to total to avoid div by zero
  sd_srad <- Rfast::rowVars(srad, std = TRUE)
  # Calculate CV * 100
  bio23V <- ifelse(mean_unit_srad <= 0, 0, (sd_srad / mean_unit_srad) * 100)
  bio23V <- cbind(bio23 = bio23V, cell = cell)
  return(bio23V)
}

#' @title bio24_fast: Radiation of Wettest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the wettest (highest precipitation).
#' @param speriod A numeric **matrix** of solar radiation values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param pperiod_max_idx An integer **vector** indicating the column index (1-based) of the wettest period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio24" (mean solar radiation of wettest period) and "cell".
#' @keywords internal
bio24_fast <- function(speriod, pperiod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(pperiod_max_idx < 1, na.rm=TRUE) || any(pperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio24_fast: Some max_prec_period indices are out of bounds.")
  }
  # Extract the mean for the wettest period
  bio24V <- speriod[cbind(seq_len(nrow(speriod)), pperiod_max_idx)]
  bio24V <- cbind(bio24 = bio24V, cell = cell)
  return(bio24V)
}

#' @title bio25_fast: Radiation of Driest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the driest (lowest precipitation).
#' @param speriod A numeric **matrix** of solar radiation values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param pperiod_min_idx An integer **vector** indicating the column index (1-based) of the driest period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio25" (mean solar radiation of driest period) and "cell".
#' @keywords internal
bio25_fast <- function(speriod, pperiod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(pperiod_min_idx < 1, na.rm=TRUE) || any(pperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio25_fast: Some min_prec_period indices are out of bounds.")
  }
  # Extract the mean for the driest period
  bio25V <- speriod[cbind(seq_len(nrow(speriod)), pperiod_min_idx)]
  bio25V <- cbind(bio25 = bio25V, cell = cell)
  return(bio25V)
}

#' @title bio26_fast: Radiation of Warmest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the warmest (highest temperature).
#' @param speriod A numeric **matrix** of solar radiation values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param tperiod_max_idx An integer **vector** indicating the column index (1-based) of the warmest period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio26" (mean solar radiation of warmest period) and "cell".
#' @keywords internal
bio26_fast <- function(speriod, tperiod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio26_fast: Some max_temp_period indices are out of bounds.")
  }
  # Extract the mean for the warmest period
  bio26V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_max_idx)]
  bio26V <- cbind(bio26 = bio26V, cell = cell)
  return(bio26V)
}

#' @title bio27_fast: Radiation of Coldest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the coldest (lowest temperature).
#' @param speriod A numeric **matrix** of solar radiation values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param tperiod_min_idx An integer **vector** indicating the column index (1-based) of the coldest period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio27" (mean solar radiation of coldest period) and "cell".
#' @keywords internal
bio27_fast <- function(speriod, tperiod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio27_fast: Some min_temp_period indices are out of bounds.")
  }
  # Extract the mean for the coldest period
  bio27V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_min_idx)]
  bio27V <- cbind(bio27 = bio27V, cell = cell)
  return(bio27V)
}

#' @title bio28_fast: Mean Moisture
#' @description Calculates mean moisture across all temporal units (usually 12 months).
#' @param mois A numeric **matrix** where **rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., 12 months).
#' @param cell A numeric or character **vector** of original cell IDs. 
#'   Its length must be exactly equal to the number of rows in `mois`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio28" (the calculated mean) and "cell" (the IDs).
#' @keywords internal
bio28_fast <- function(mois, cell) {
  bio28V <- Rfast::rowmeans(mois)
  bio28V <- cbind(bio28 = bio28V, cell = cell)
  return(bio28V)
}

#' @title bio29_fast: Highest Moisture Unit
#' @description Identifies the highest moisture of the temporal unit with the highest value. 
#'   If `index_vector` is `NULL`, it calculates the row-wise maximum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @param mois A numeric **matrix** of moisture values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `mois`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `mois`. 
#'   Values must be between 1 and `ncol(mois)`. 
#'   This is typically used to extract the moisture of a specific unit identified by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio29" (the maximum moisture) and "cell".
#' @keywords internal
bio29_fast <- function(mois, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(mois)) stop("bio29_fast: Length mismatch: index_vector vs mois rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(mois)
     valid_rows <- which(!is_invalid_idx)
     bio29V <- rep(NA_real_, nrow(mois))
     if(length(valid_rows) > 0) {
        bio29V[valid_rows] <- mois[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio29_fast: Some static indices were NA or out of bounds.")
  } else {
    bio29V <- Rfast::rowMaxs(mois, value = TRUE)
  }
  bio29V <- cbind(bio29 = bio29V, cell = cell)
  return(bio29V)
}

#' @title bio30_fast: Lowest Moisture Unit
#' @description Identifies the lowest moisture of the temporal unit with the lowest value. 
#'   If `index_vector` is `NULL`, it calculates the row-wise minimum. 
#'   If `index_vector` is provided, it extracts the value from the specific column index for each row.
#' @param mois A numeric **matrix** of moisture values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units (e.g., months).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `mois`.
#' @param index_vector (Optional) An integer **vector** of column indices (1-based). 
#'   If provided, its length must be exactly equal to the number of rows in `mois`. 
#'   Values must be between 1 and `ncol(mois)`. 
#'   This is typically used to extract the moisture of a specific unit identified by another metric.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio30" (the minimum moisture) and "cell".
#' @keywords internal
bio30_fast <- function(mois, cell, index_vector = NULL) {
  if (!is.null(index_vector)) {
     if (length(index_vector) != nrow(mois)) stop("bio30_fast: Length mismatch: index_vector vs mois rows.")
     is_invalid_idx <- is.na(index_vector) | index_vector < 1 | index_vector > ncol(mois)
     valid_rows <- which(!is_invalid_idx)
     bio30V <- rep(NA_real_, nrow(mois))
     if(length(valid_rows) > 0) {
        bio30V[valid_rows] <- mois[cbind(valid_rows, index_vector[valid_rows])]
     }
     if(any(is_invalid_idx)) warning("bio30_fast: Some static indices were NA or out of bounds.")
  } else {
    bio30V <- Rfast::rowMins(mois, value = TRUE)
  }
  bio30V <- cbind(bio30 = bio30V, cell = cell)
  return(bio30V)
}

#' @title bio31_fast: Moisture Seasonality (Std Dev * 100)
#' @description Calculates Moisture Seasonality, defined as the standard deviation of moisture values across all temporal units (e.g., 12 months), multiplied by 100.
#' @param mois A numeric **matrix** containing moisture values. **Rows** represent spatial units (cells) 
#'   and **columns** represent temporal units.
#' @param n_units A single **integer** representing the number of temporal units (e.g., 12).
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `mois`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells. 
#'   The columns are named "bio31" (Moisture Seasonality) and "cell".
#' @keywords internal
bio31_fast <- function(mois, n_units, cell) {
  bio31V <- Rfast::rowVars(mois, std = TRUE) * 100
  bio31V <- cbind(bio31 = bio31V, cell = cell)
  return(bio31V)
}

#' @title bio32_fast: Mean Moisture of Most Moist Period
#' @description Calculates the mean moisture of the specific rolling period identified as the most moist (highest moisture).
#' @param speriod A numeric **matrix** of moisture values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param speriod_max_idx An integer **vector** indicating the column index (1-based) of the most moist period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio32" (mean moisture of most moist period) and "cell".
#' @keywords internal
bio32_fast <- function(speriod, speriod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(speriod_max_idx < 1, na.rm=TRUE) || any(speriod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio32_fast: Some max_prec_period indices are out of bounds.")
  }
  # Extract the mean for the most moist period
  bio32V <- speriod[cbind(seq_len(nrow(speriod)), speriod_max_idx)]
  bio32V <- cbind(bio32 = bio32V, cell = cell)
  return(bio32V)
}

#' @title bio33_fast: Mean Moisture of Least Moist Period
#' @description Calculates the mean moisture of the specific rolling period identified as the least moist (lowest moisture).
#' @param speriod A numeric **matrix** of moisture values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param speriod_min_idx An integer **vector** indicating the column index (1-based) of the least moist period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio33" (mean moisture of least moist period) and "cell".
#' @keywords internal
bio33_fast <- function(speriod, speriod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(speriod_min_idx < 1, na.rm=TRUE) || any(speriod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio33_fast: Some min_prec_period indices are out of bounds.")
  }
  # Extract the mean for the least moist period
  bio33V <- speriod[cbind(seq_len(nrow(speriod)), speriod_min_idx)]
  bio33V <- cbind(bio33 = bio33V, cell = cell)
  return(bio33V)
}

#' @title bio34_fast: Mean Moisture of Warmest Period
#' @description Calculates the mean moisture of the specific rolling period identified as the warmest (highest temperature).
#' @param speriod A numeric **matrix** of moisture values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param tperiod_max_idx An integer **vector** indicating the column index (1-based) of the warmest period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio34" (mean moisture of warmest period) and "cell".
#' @keywords internal
bio34_fast <- function(speriod, tperiod_max_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_max_idx < 1, na.rm=TRUE) || any(tperiod_max_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio34_fast: Some max_temp_period indices are out of bounds.")
  }
  # Extract the mean for the warmest period
  bio34V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_max_idx)]
  bio34V <- cbind(bio34 = bio34V, cell = cell)
  return(bio34V)
}

#' @title bio35_fast: Mean Moisture of Coldest Period
#' @description Calculates the mean moisture of the specific rolling period identified as the coldest (lowest temperature).
#' @param speriod A numeric **matrix** of moisture values (means) for each rolling period. 
#'   **Rows** represent spatial units (cells). **Columns** represent the rolling periods. 
#' @param tperiod_min_idx An integer **vector** indicating the column index (1-based) of the coldest period for each row. 
#'   Its length must be exactly equal to the number of rows in `speriod`.
#' @param cell A vector of original cell IDs. Its length must be exactly equal to the number of rows in `speriod`.
#' @return A **matrix** with dimensions `c(N, 2)`, where N is the number of input cells.
#'   The columns are named "bio35" (mean moisture of coldest period) and "cell".
#' @keywords internal
bio35_fast <- function(speriod, tperiod_min_idx, cell) {
  num_period_cols <- ncol(speriod) - 2
  if (any(tperiod_min_idx < 1, na.rm=TRUE) || any(tperiod_min_idx > num_period_cols, na.rm=TRUE)) {
    warning("bio35_fast: Some min_temp_period indices are out of bounds.")
  }
  # Extract the mean for the coldest period
  bio35V <- speriod[cbind(seq_len(nrow(speriod)), tperiod_min_idx)]
  bio35V <- cbind(bio35 = bio35V, cell = cell)
  return(bio35V)
}