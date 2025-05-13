# Helper functions
# -----------------------

# Average Temperature
#' @export
t_avg <- function(tmin, tmax) {
  misqua <- (tmin + tmax) / 2
  names(misqua) <- paste0("tavg_", 1:terra::nlyr(misqua))
  return(misqua)
}

# CV terra
#' @export
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
#' @export
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

#' @export
testGeom <- function(x, y) {
  testGeom <- terra::compareGeom(x, y, lyrs = TRUE)
  if (testGeom == TRUE) {
    return(x)
  }
}

#' @export
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
