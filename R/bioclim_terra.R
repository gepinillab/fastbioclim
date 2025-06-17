#' In-Memory Bioclimatic Variable Calculation
#'
#' Internal function to calculate bioclimatic variables using `terra` functions.
#' It is designed for datasets that can fit into RAM.
#'
#' @param bios Numeric vector of variables to compute.
#' @param ... `SpatRaster` objects for climate variables (e.g., `tmin`, `tmax`)
#'   and static indices (e.g., `warmest_period`).
#' @param output_dir Character, path to save final rasters.
#' @param period_length Integer, length of a calculation period.
#' @param circular Logical, whether to wrap periods.
#' @param gdal_opt Character vector of GDAL options for writing.
#' @param overwrite Logical, whether to overwrite existing files.
#' @return A `terra::SpatRaster` object pointing to the newly created files.
#' @keywords internal
#' @seealso The user-facing wrapper function `derive_bioclim()`.
#' 
bioclim_terra <- function(bios, 
  tmin = NULL, 
  tmax = NULL, 
  tavg = NULL, 
  prcp = NULL,
  srad = NULL, 
  mois = NULL, 
  period_length = 3,
  circular = TRUE, 
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE, 
  output_dir = tempdir(),
  ...) {
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
  }
  
  # Bios that requires solar radiation
  req_srad <- c(20, 21, 22, 23, 24, 25, 26, 27)
  if (any(req_srad %in% bios)) {
    if (is.null(srad)) {
      stop(paste0(paste0("Bio", sprintf("%02d", req_srad[req_srad %in% bios]),
                         collapse = ", "),
                  " require(s) srad."))
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
    wet <- fastbioclim::get_window(prcp, period_length, circular)
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
    tmp <- fastbioclim::get_window(tavg, period_length, circular) / period_length
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
    prad <- fastbioclim::get_window(srad, period_length, circular) / period_length
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
    pmois <- fastbioclim::get_window(mois, period_length, circular) / period_length
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
             " was(were) built with a period of ", period_length,
             " units with", if(circular == FALSE) "out", " circularity.")
    )
  }

  # Create a unique spatRaster
  bios_rast <- terra::rast(mget(paste0("bio", sprintf("%02d", bios))))
  output_files <- file.path(output_dir, paste0(names(bios_rast), ".tif"))

  message("Writing GeoTIFFs...")
  terra::writeRaster(
    bios_rast,
    filename = output_files,
    overwrite = overwrite,
    gdal = gdal_opt
  )
  return(output_files)
}
