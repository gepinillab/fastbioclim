#' @title Get several bioclimatic variables
#' @description
#' The clima function generates bioclimatic variables commonly used in ecological 
#' and environmental analyses. This function builds upon the biovars function from 
#' the dismo package and adheres to the methodologies outlined in the ANUCLIM 6.1 manual 
#' and O’Donnell & Ignizio (2012). It calculates a comprehensive set of 35 
#' bioclimatic variables based on input rasters for minimum temperature, 
#' maximum temperature, average temperature, precipitation, solar radiation, 
#' and moisture. Users can specify which bioclimatic variables to compute 
#' using the bios parameter. The function supports customizable temporal 
#' aggregation (e.g., quarters or semesters) and includes options for circular 
#' period calculations. Input consistency is validated to ensure matching 
#' extents, resolutions, and NA handling, making it robust for large-scale 
#' spatial applications. Outputs are returned as a SpatRaster object 
#' containing the selected bioclimatic layers, ready for downstream analysis.
#' 
#' @param bios Numeric vector with bios number
#' @param tmin spatRaster with minimum temperature
#' @param tmax spatRaster with maximum temperature
#' @param tavg spatRaster with average temperature. It could be NULL if tmin and
#'  tmax are provided
#' @param prcp spatRaster with precipitation
#' @param srad spatRaster with solar radiation
#' @param mois spatRaster with moisture
#' @param period numeric. Length of period to summarize data (e.g., quarters,
#' semesters). If using monthly data, a quarter (3-months) will be used to
#' calculate Bios 08, 09, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 33, 34, 35.
#' @param circular TRUE/FALSE. Calculate periods that include the first and last
#' units. For example, if using mean monthly data and quarters, `circular=TRUE` will also
#' calculate Nov-Dec-Jan and Dec-Jan-Feb.
#' @param checkNA TRUE/FALSE. Check for unexpected NA values in the input data.
#' @param stopNA TRUE/FALSE. Stop execution if unexpected NA values are detected.
#' @param ... Additional arguments to define spatRasters that define a static 
#' period or unit index for variables calculation could be useful to keep same 
#' index along a time series or define a static season of interest, including:
#'   \describe{
#'     \item{warmest_period}{Description of param1 (default: NULL)}
#'     \item{coldest_period}{Description of param1 (default: NULL)}
#'     \item{wettest_period}{Description of param1 (default: NULL)}
#'     \item{driest_period}{Description of param1 (default: NULL)}
#'     \item{high_mois_period}{Description of param1 (default: NULL)}
#'     \item{low_mois_period}{Description of param1 (default: NULL)}
#'     \item{warmest_unit}{Description of param1 (default: NULL)}
#'     \item{coldest_unit}{Description of param1 (default: NULL)}
#'     \item{wettest_unit}{Description of param1 (default: NULL)}
#'     \item{driest_unit}{Description of param1 (default: NULL)}
#'     \item{high_rad_unit}{Description of param1 (default: NULL)}
#'     \item{low_rad_unit}{Description of param1 (default: NULL)}
#'     \item{high_mois_unit}{Description of param1 (default: NULL)}
#'     \item{low_mois_unit}{Description of param1 (default: NULL)}
#'   }
#' 
#' @return An SpatRaster with 35 bioclimatic variables or a subset of them:
#' \describe{
#'   \item{bio01}{Mean Temperature of Units}
#'   \item{bio02}{Mean Diurnal Range}
#'   \item{bio03}{Isothermality}
#'   \item{bio04}{Temperature Seasonality}
#'   \item{bio05}{Max Temperature of Warmest Unit}
#'   \item{bio06}{Min Temperature of Coldest Unit}
#'   \item{bio07}{Temperature Range of Units}
#'   \item{bio08}{Mean Temperature of Wettest Period}
#'   \item{bio09}{Mean Temperature of Driest Period}
#'   \item{bio10}{Mean Temperature of Warmest Period}
#'   \item{bio11}{Mean Temperature of Coldest Period}
#'   \item{bio12}{Precipitation Sum}
#'   \item{bio13}{Precipitation of Wettest Unit}
#'   \item{bio14}{Precipitation of Driest Unit}
#'   \item{bio15}{Precipitation Seasonality}
#'   \item{bio16}{Precipitation of Wettest Period}
#'   \item{bio17}{Precipitation of Driest Period}
#'   \item{bio18}{Precipitation of Warmest Period}
#'   \item{bio19}{Precipitation of Coldest Period}
#'   \item{bio20}{Mean Radiation of Units}
#'   \item{bio21}{Highest Radiation Unit}
#'   \item{bio22}{Lowest Radiation Unit}
#'   \item{bio23}{Radiation Seasonality}
#'   \item{bio24}{Radiation of Wettest Period}
#'   \item{bio25}{Radiation of Driest Period}
#'   \item{bio26}{Radiation of Warmest Period}
#'   \item{bio27}{Radiation of Coldest Period}
#'   \item{bio28*}{Mean Moisture Content Of Units}
#'   \item{bio29*}{Highest Moisture Content Unit}
#'   \item{bio30*}{Lowest Moisture Content Unit}
#'   \item{bio31*}{Moisture Content Seasonality}
#'   \item{bio32*}{Mean Moisture Content of Most Moist Period}
#'   \item{bio33*}{Mean Moisture Content of Least Moist Period}
#'   \item{bio34*}{Mean Moisture Content of Warmest Period}
#'   \item{bio35*}{Mean Moisture Content of Coldest Period}
#' }
#' 
#' @note 
#' *The original moisture variables proposed in the ANUCLIM manual are based 
#' on the Moisture Index (MI). However, this function allows users to calculate 
#' moisture-based bioclimatic variables using other units of moisture 
#' as inputs, offering greater flexibility in input data usage.
#' 
#' @references
#' O’Donnell, M. S., & Ignizio, D. A. (2012). Bioclimatic predictors for supporting ecological applications in the conterminous United States (Vol. 691).  
#' ANUCLIM 6.1 User Guide. Centre for Resource and Environmental Studies, The Australian National University.  
#' Hijmans, R. J., Phillips, S., Leathwick, J., & Elith, J. (2017). `dismo`: Species Distribution Modeling. R package version 1.1-4.
#' 
#' @export
#'
clima <- function(bios, tmin = NULL, tmax = NULL, tavg = NULL, prcp = NULL,
                  srad = NULL, mois = NULL, period = 3, circular = TRUE, 
                  checkNA = TRUE, stopNA = TRUE, ...) {
  # Add dor arguments
  dot_args <- list(...)
  if ("warmest_unit" %in% names(dot_args)) warmest_unit <- dot_args$warmest_unit
  if ("coldest_unit" %in% names(dot_args)) coldest_unit <- dot_args$coldest_unit
  if ("wettest_unit" %in% names(dot_args)) wettest_unit <- dot_args$wettest_unit
  if ("driest_unit" %in% names(dot_args)) driest_unit <- dot_args$driest_unit
  if ("high_rad_unit" %in% names(dot_args)) high_rad_unit <- dot_args$high_rad_unit
  if ("low_rad_unit" %in% names(dot_args)) low_rad_unit <- dot_args$low_rad_unit
  if ("high_mois_unit" %in% names(dot_args)) high_mois_unit <- dot_args$high_mois_unit
  if ("low_mois_unit" %in% names(dot_args)) low_mois_unit <- dot_args$low_mois_unit
  if ("warmest_period" %in% names(dot_args)) warmest_period <- dot_args$warmest_period
  if ("coldest_period" %in% names(dot_args)) coldest_period <- dot_args$coldest_period
  if ("wettest_period" %in% names(dot_args)) wettest_period <- dot_args$wettest_period
  if ("driest_period" %in% names(dot_args)) driest_period <- dot_args$driest_period
  if ("high_mois_period" %in% names(dot_args)) high_mois_period <- dot_args$high_mois_period
  if ("low_mois_period" %in% names(dot_args)) low_mois_period <- dot_args$low_mois_period
  # Check for same extent, number of rows and columns, projection,
  # resolution, and origin
  sameGeom <- class(purrr::reduce(list(tmin, tmax, tavg, prcp, srad, mois) |>
                                    purrr::discard(is.null),
                                  fastbioclim::testGeom))
  if (sameGeom == "SpatRaster") {
    message("SpatRasters have same extent, number of rows and columns, ",
            "projection, resolution, and origin")
  }

  # CHECKING INPUTS AVAILABLE
  # Bios that requires prcp
  req_prcp <- c(8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 24, 25)
  if (exists("wettest_period")) {
    if (8 %in% req_prcp) req_prcp <- req_prcp[which(req_prcp != 8)]
    if (24 %in% req_prcp) req_prcp <- req_prcp[which(req_prcp != 24)] 
  }
  if (exists("driest_period")) {
    if (9 %in% req_prcp) req_prcp <- req_prcp[which(req_prcp != 9)]
    if (25 %in% req_prcp) req_prcp <- req_prcp[which(req_prcp != 25)] 
  }
  if (any(req_prcp %in% bios)) {
    if (is.null(prcp)) {
      stop(paste0(paste0("Bio", sprintf("%02d", req_prcp[req_prcp %in% bios]),
                         collapse = ", "),
                  " require(s) prcp."))
    }
    if (checkNA == TRUE) {
      prcp_na <- mismatch_NA(prcp)
      if (prcp_na$logical == TRUE & stopNA == TRUE) {
        stop("prcp has unexpected NA values")
      }
      prcp_sum <- prcp_na$sum_lyr
    }
  }

  # Bios that requires temperature
  req_temp <- c(1, 2, 3, 4, 7, 8, 9, 10, 11, 18, 19, 26, 27, 34, 35)
  if (exists("warmest_period")) {
    if (18 %in% req_temp) req_temp <- req_temp[which(req_temp != 18)]
    if (26 %in% req_temp) req_temp <- req_temp[which(req_temp != 26)]
    if (34 %in% req_temp) req_temp <- req_temp[which(req_temp != 34)]
  }
  if (exists("coldest_period")) {
    if (19 %in% req_temp) req_temp <- req_temp[which(req_temp != 19)]
    if (27 %in% req_temp) req_temp <- req_temp[which(req_temp != 27)]
    if (35 %in% req_temp) req_temp <- req_temp[which(req_temp != 35)]
  }
  if (any(c(5, 6, req_temp) %in% bios)) {
    # Bios that requires just tmax
    if (5 %in% bios & is.null(tmax)) stop("Bio05 requires tmax.")
    # Bios that requires just tmin
    if (6 %in% bios & is.null(tmin)) stop("Bio06 requires tmin.")
    # Bios that requires tmin and tmax
    if (is.null(tmin) & is.null(tmax)) {
      stop(paste0(paste0("Bio", sprintf("%02d", req_temp[req_temp %in% bios]),
                         collapse = ", "),
                  " require(s) tmin and tmax."))
    }
    # Check for NA (tmax)
    if (any(c(5, req_temp) %in% bios) & checkNA == TRUE) {
      tmax_na <- mismatch_NA(tmax)
      if (tmax_na$logical == TRUE & stopNA == TRUE) {
        stop("tmax has unexpected NA values")
      }
      tmax_sum <- tmax_na$sum_lyr
    }
    # Check for NA (tmin)
    if (any(c(6, req_temp) %in% bios) & checkNA == TRUE) {
      tmin_na <- mismatch_NA(tmin)
      if (tmin_na$logical == TRUE & stopNA == TRUE) {
        stop("tmin has unexpected NA values")
      }
      tmin_sum <- tmin_na$sum_lyr
    }
  }
  
  # Bios that requires solar radiation
  req_srad <- c(20, 21, 22, 23, 24, 25, 26, 27)
  if (any(req_srad %in% bios)) {
    if (is.null(srad)) {
      stop(paste0(paste0("Bio", sprintf("%02d", req_srad[req_srad %in% bios]),
                         collapse = ", "),
                  " require(s) srad."))
    }
    if (checkNA == TRUE) {
      srad_na <- mismatch_NA(srad)
      if (srad_na$logical == TRUE & stopNA == TRUE) {
        stop("srad has unexpected NA values")
      }
      srad_sum <- srad_na$sum_lyr
    }
  }
  
  # Bios that requires moisture
  req_mois <- c(28, 29, 30, 31, 32, 33, 34, 35)
  if (any(req_mois %in% bios)) {
    if (is.null(mois)) {
      stop(paste0(paste0("Bio", sprintf("%02d", req_mois[req_mois %in% bios]),
                         collapse = ", "),
                  " require(s) mois."))
    }
    if (checkNA == TRUE) {
      mois_na <- mismatch_NA(mois)
      if (mois_na$logical == TRUE & stopNA == TRUE) {
        stop("mois has unexpected NA values")
      }
      mois_sum <- mois_na$sum_lyr
    }
  }

  if (checkNA == TRUE) {
    # Sum all rasters
    if (!exists("tmin_sum")) tmin_sum <-  NULL
    if (!exists("tmax_sum")) tmax_sum <-  NULL
    if (!exists("prcp_sum")) prcp_sum <-  NULL
    if (!exists("srad_sum")) srad_sum <-  NULL
    if (!exists("mois_sum")) mois_sum <-  NULL

    intra_na <- purrr::reduce(list(tmin_sum, tmax_sum, prcp_sum, srad_sum, mois_sum) |>
                                 purrr::discard(is.null), c) |>
      sum() |>
      terra::unique() |>
      unlist()
    # Delete values that shared pixel NA in all layers
    intra_na <- intra_na[intra_na != 0]
    # Check if there is sum of pixel with values is equal
    if (length(intra_na) != 1) {
      if (stopNA == TRUE) {
        stop("SpatRaster don't share same NA values")
      } else {
        message("SpatRasters (tmin, tmax, prcp, srad and/or mois) don't share same NAs values")
      }
    }
  }

  # Bios that requires tavg
  if (any(c(1, 4, 8, 9, 10, 11, 18, 19) %in% bios)) {
    if (is.null(tavg)) {
      if (is.null(tmin) | is.null(tmax)) {
        stop("tavg cannot be calculated because tmin and/or tmax are NULL")
      } else {
        tavg <- fastbioclim::t_avg(tmin, tmax)
      }
    } 
  }

  ## ONLY TEMPERATURE
  # Bio01
  if (1 %in% bios) bio01 <- bio01_terra(tavg)
  # Bio02
  if (any(2:3 %in% bios)) bio02 <- bio02_terra(tmin, tmax)
  # Bio04
  if (4 %in% bios) bio04 <- bio04_terra(tavg)
  # Bio05
  if (any(c(3, 5, 7) %in% bios)) bio05 <- bio05_terra(tmax, ...)
  # Bio06
  if (any(c(3, 6, 7) %in% bios)) bio06 <- bio06_terra(tmin, ...)
  # Bio07
  if (any(c(3, 7) %in% bios)) bio07 <- bio07_terra(bio05, bio06)
  # Bio03
  if (3 %in% bios) bio03 <- bio03_terra(bio02, bio07)
  
  ## ONLY PRECIPITATION
  # Bio12
  if (12 %in% bios) bio12 <- bio12_terra(prcp)
  # Bio13
  if (13 %in% bios) bio13 <- bio13_terra(prcp, ...)
  # Bio14
  if (14 %in% bios) bio14 <- bio14_terra(prcp, ...)
  # Bio15
  if (15 %in% bios) bio15 <- bio15_terra(prcp)
  
  
  ## ONLY PRECIPITATION PERIOD
  if (any(c(8:9, 16:19, 24:25) %in% bios)) {
    wet <- fastbioclim::get_window(prcp, period, circular)
    if (any(c(8, 16, 24) %in% bios) & !exists("wettest_period")) {
      wettest_period <- terra::which.max(wet)
    }
    if (any(c(9, 17, 25) %in% bios) & !exists("driest_period")) {
      driest_period <- terra::which.min(wet)
    }
  }
  # Bio16
  if (16 %in% bios) bio16 <- bio16_terra(wet, wettest_period)
  # Bio17
  if (17 %in% bios) bio17 <- bio17_terra(wet, driest_period)
  
  ### ONLY TEMPERATURE PERIOD
  
  if (any(c(8:11, 18:19, 26:27, 34:35) %in% bios)) {
    if (any(c(8:9) %in% bios) | 
        (any(c(10, 18, 26, 34) %in% bios) & !exists("warmest_period")) |
        (any(c(11, 19, 27, 35) %in% bios) & !exists("coldest_period"))) {
      tmp <- fastbioclim::get_window(tavg, period, circular) / period
    }
    if (any(c(10, 18, 26, 34) %in% bios) & !exists("warmest_period")) {
      warmest_period <- terra::which.max(tmp)
    }
    if (any(c(11, 19, 27, 35) %in% bios) & !exists("coldest_period")) {
      coldest_period <- terra::which.min(tmp)
    }
  }
  # Bio10
  if (10 %in% bios) bio10 <- bio10_terra(tmp, warmest_period)
  # Bio11
  if (11 %in% bios) bio11 <- bio11_terra(tmp, coldest_period)
  
  ## ONLY SOLAR RADIATION
  # Bio20
  if (20 %in% bios) bio20 <- bio20_terra(srad)
  # Bio21
  if (21 %in% bios) bio21 <- bio21_terra(srad, ...)
  # Bio22
  if (22 %in% bios) bio22 <- bio22_terra(srad, ...)
  # Bio23
  if (23 %in% bios) bio23 <- bio23_terra(srad)
  
  ### GET SOLAR RADIATION PERIOD
  if (any(c(24:27) %in% bios)) {
    prad <- fastbioclim::get_window(srad, period, circular) / period
  }
  
  ## ONLY MOISTURE
  # Bio28
  if (28 %in% bios) bio28 <- bio28_terra(mois)
  # Bio29
  if (29 %in% bios) bio29 <- bio29_terra(mois, ...)
  # Bio30
  if (30 %in% bios) bio30 <- bio30_terra(mois, ...)
  # Bio31
  if (31 %in% bios) bio31 <- bio31_terra(mois)
  
  ### ONLY MOISTURE PERIOD
  if (any(c(32:35) %in% bios)) {
    pmois <- fastbioclim::get_window(mois, period, circular) / period
    if ((32 %in% bios) & !exists("high_mois_period")) {
      high_mois_period <- terra::which.max(pmois)
    }
    if ((33 %in% bios) & !exists("low_mois_period")) {
      low_mois_period <- terra::which.min(pmois)
    }
  }
  # Bio32
  if (32 %in% bios) bio32 <- bio32_terra(pmois, high_mois_period)
  # Bio33
  if (33 %in% bios) bio33 <- bio33_terra(pmois, low_mois_period)
  
  #### COMBINED PERIODS
  # Bio08
  if (8 %in% bios) bio08 <- bio08_terra(tmp, wettest_period)
  # Bio09
  if (9 %in% bios) bio09 <- bio09_terra(tmp, driest_period)
  # Bio18
  if (18 %in% bios) bio18 <- bio18_terra(wet, warmest_period)
  # Bio19
  if (19 %in% bios) bio19 <- bio19_terra(wet, coldest_period)
  # Bio24
  if (24 %in% bios) bio24 <- bio24_terra(prad, wettest_period)
  # Bio25
  if (25 %in% bios) bio25 <- bio25_terra(prad, driest_period)
  # Bio26
  if (26 %in% bios) bio26 <- bio26_terra(prad, warmest_period)
  # Bio27
  if (27 %in% bios) bio27 <- bio27_terra(prad, coldest_period)
  # Bio34
  if (34 %in% bios) bio34 <- bio34_terra(pmois, warmest_period)
  # Bio35
  if (35 %in% bios) bio35 <- bio35_terra(pmois, coldest_period)
  
  # Window message
  if (any(c(8:11, 16:19, 24:27, 32:35) %in% bios)) {
    message(
      paste0(paste0("Bio", sprintf("%02d", c(8:11, 16:19, 24:27, 32:35)[c(8:11, 16:19, 24:27, 32:35) %in% bios]),
                    collapse = ", "),
             " was(were) built with a period of ", period,
             " units with", if(circular == FALSE) "out", " circularity.")
    )
  }

  # Create a unique spatRaster
  bios_rast <- terra::rast(mget(paste0("bio", sprintf("%02d", bios))))
  return(bios_rast)
}
