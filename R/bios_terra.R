# Function to create analogous bioclimatic variables
# ----------
#' @title bio01_terra: Mean Temperature of Units
#' @description Calculates mean temperature across all temporal units.
#' @param tavg spatRaster of average temperatures for each unit.
#' @return spatRaster with "bio01".
#' @keywords internal
bio01_terra <- function(tavg) {
  ata <- terra::app(tavg, mean, na.rm = TRUE)
  names(ata) <- "bio01"
  return(ata)
}

#' @title bio02_terra: Mean Diurnal Range
#' @description Calculates the mean of (tmax - tmin) across all temporal units.
#' @param tmin spatRaster of minimum temperatures for each unit.
#' @param tmax spatRaster of maximum temperatures for each unit.
#' @return spatRaster with "bio02".
#' @keywords internal
bio02_terra <- function(tmin, tmax) {
  bosa <- terra::app(tmax - tmin, mean, na.rm = TRUE)
  names(bosa) <- "bio02"
  return(bosa)
}

#' @title bio03_terra: Isothermality
#' @description Calculates (bio02 / bio07) * 100.
#' @param bio02 spatRaster of bio02 values.
#' @param bio07 spatRaster of bio07 values.
#' @return spatRaster with "bio03".
#' @keywords internal
bio03_terra <- function(bio02, bio07) {
   mica <- 100 * bio02 / bio07
  names(mica) <- "bio03"
  return(mica)
}

#' @title bio04_terra: Temperature Seasonality (Std Dev * 100)
#' @description Calculates the standard deviation of average temperatures across units, multiplied by 100.
#' @param tavg spatRaster of average temperatures for each unit.
#' @return spatRaster with "bio04".
#' @keywords internal
bio04_terra <- function(tavg) {
  muihica <- 100 * terra::stdev(tavg, pop = FALSE, na.rm = TRUE)
  names(muihica) <- "bio04"
  return(muihica)
}

#' @title bio05_terra: Max Temperature of Warmest Unit
#' @description Identifies max temperature of the warmest unit, potentially using a static index.
#' @param tmax spatRaster of maximum temperatures for each unit.
#' @return spatRaster with "bio05".
#' @keywords internal
bio05_terra <- function(tmax, warmest_unit = NULL) {
  if (!is.null(warmest_unit)) {
    hisca <- terra::selectRange(tmax, warmest_unit)
  } else {
    hisca <- terra::app(tmax, max, na.rm = TRUE)
  }
  names(hisca) <- "bio05"
  return(hisca)
}

#' @title bio06_terra: Min Temperature of Coldest Unit
#' @description Identifies min temperature of the coldest unit, potentially using a static index.
#' @param tmin spatRaster of minimum temperatures for each unit.
#' @return spatRaster with "bio06".
#' @keywords internal
bio06_terra <- function(tmin, coldest_unit = NULL) {
  if (!is.null(coldest_unit)) {
    ta <- terra::selectRange(tmin, coldest_unit)
  } else {
    ta <- terra::app(tmin, min, na.rm = TRUE)
  }
  names(ta) <- "bio06"
  return(ta)
}

#' @title bio07_terra: Temperature Annual Range (bio05 - bio06)
#' @description Calculates the difference between bio05 and bio06.
#' @param bio05 spatRaster of bio05 values.
#' @param bio06 spatRaster of bio06 values.
#' @return spatRaster with "bio07".
#' @keywords internal
bio07_terra <- function(bio05, bio06) {
  cuhupcua <- bio05 - bio06
  names(cuhupcua) <- "bio07"
  return(cuhupcua)
}

#' @title bio08_terra: Mean Temperature of Wettest Period
#' @description Calculates mean temperature of the period with the highest precipitation sum.
#' @param tmp spatRaster of temperature period sums.
#' @param wettest_period spatRaster indicating the index (1-based) of the wettest period for each cell.
#' @return spatRaster with "bio08".
#' @keywords internal
bio08_terra <- function(tmp, wettest_period) {
  suhusa <- terra::selectRange(tmp, wettest_period)
  names(suhusa) <- "bio08"
  return(suhusa)
}

#' @title bio09_terra: Mean Temperature of Driest Period
#' @description Calculates mean temperature of the period with the lowest precipitation sum.
#' @param tmp spatRaster of temperature period sums.
#' @param driest_period Vector indicating the index (1-based) of the driest period.
#' @return spatRaster with "bio09".
#' @keywords internal
bio09_terra <- function(tmp, driest_period) {
  aca <- terra::selectRange(tmp, driest_period)
  names(aca) <- "bio09"
  return(aca)
}

#' @title bio10_terra: Mean Temperature of Warmest Period
#' @description Calculates mean temperature of the period with the highest temperature sum.
#' @param tmp spatRaster of temperature period sums.
#' @param warmest_period Vector indicating the index (1-based) of the warmest period.
#' @return spatRaster with "bio10".
#' @keywords internal
bio10_terra <- function(tmp, warmest_period) {
  ubchihica <- terra::selectRange(tmp, warmest_period)
  names(ubchihica) <- "bio10"
  return(ubchihica)
}

#' @title bio11_terra: Mean Temperature of Coldest Period
#' @description Calculates mean temperature of the period with the lowest temperature sum.
#' @param tmp spatRaster of temperature period sums.
#' @param coldest_period Vector indicating the index (1-based) of the coldest period.
#' @return spatRaster with "bio11".
#' @keywords internal
bio11_terra <- function(tmp, coldest_period) {
  quihicha_ata <- terra::selectRange(tmp, coldest_period)
  names(quihicha_ata) <- "bio11"
  return(quihicha_ata)
}

#' @title bio12_terra: Total Precipitation
#' @description Calculates the sum of precipitation values across all units.
#' @param prcp spatRaster of precipitation values for each unit.
#' @return spatRaster with "bio12".
#' @keywords internal
bio12_terra <- function(prcp) {
  quihicha_bosa <- terra::app(prcp, sum, na.rm = TRUE)
  names(quihicha_bosa) <- "bio12"
  return(quihicha_bosa)
}

#' @title bio13_terra: Precipitation of Wettest Unit
#' @description Identifies precipitation of the wettest unit, potentially using a static index.
#' @param prcp spatRaster of precipitation values for each unit.
#' @return spatRaster with "bio13".
#' @keywords internal
bio13_terra <- function(prcp, wettest_unit = NULL) {
  if (!is.null(wettest_unit)) {
    quihicha_mica <- terra::selectRange(prcp, wettest_unit)
  } else {
    quihicha_mica <- terra::app(prcp, max, na.rm = TRUE)
  }
  names(quihicha_mica) <- "bio13"
  return(quihicha_mica)
}

#' @title bio14_terra: Precipitation of Driest Unit
#' @description Identifies precipitation of the driest unit, potentially using a static index.
#' @param prcp spatRaster of precipitation values for each unit.
#' @return spatRaster with "bio14".
#' @keywords internal
bio14_terra <- function(prcp, driest_unit = NULL) {
  if (!is.null(driest_unit)) {
    quihicha_muihica <- terra::selectRange(prcp, driest_unit)
  } else {
    quihicha_muihica <- terra::app(prcp, min, na.rm = TRUE)
  }
  names(quihicha_muihica) <- "bio14"
  return(quihicha_muihica)
}
#' @title bio15_terra: Precipitation Seasonality (CV)
#' @description Calculates coefficient of variation in precipitation across units.
#' @param prcp Matrix containing precipitation values for each unit.
#' @note The "1 +" is to avoid strange CVs for areas where mean rainfaill is < 1)
#' @return spatRaster with "bio15".
#' @keywords internal
bio15_terra <- function(prcp) {
  quihicha_hisca <- cv_cli(prcp)
  names(quihicha_hisca) <- "bio15"
  return(quihicha_hisca)
}

#' @title bio16_terra: Precipitation of Wettest Period
#' @description Calculates precipitation sum of the period with the highest precipitation sum.
#' @param wet spatRaster of precipitation period sums.
#' @param wettest_period Vector indicating the index (1-based) of the wettest period.
#' @return spatRaster with "bio16".
#' @keywords internal
bio16_terra <- function(wet, wettest_period) {
  quihicha_ta <- terra::selectRange(wet, wettest_period)
  names(quihicha_ta) <- "bio16"
  return(quihicha_ta)
}

#' @title bio17_terra: Precipitation of Driest Period
#' @description Calculates precipitation sum of the period with the lowest precipitation sum.
#' @param wet spatRaster of precipitation period sums.
#' @param driest_period Vector indicating the index (1-based) of the driest period.
#' @return spatRaster with "bio17".
#' @keywords internal
bio17_terra <- function(wet, driest_period) {
  quihicha_cuhupcua <- terra::selectRange(wet, driest_period)
  names(quihicha_cuhupcua) <- "bio17"
  return(quihicha_cuhupcua)
}

#' @title bio18_terra: Precipitation of Warmest Period
#' @description Calculates precipitation sum of the period with the highest temperature sum.
#' @param wet spatRaster of precipitation period sums.
#' @param warmest_period Vector indicating the index (1-based) of the warmest period.
#' @return spatRaster with "bio18".
#' @keywords internal
bio18_terra <- function(wet, warmest_period) {
  quihicha_suhusa <- terra::selectRange(wet, warmest_period)
  names(quihicha_suhusa) <- "bio18"
  return(quihicha_suhusa)
}

#' @title bio19_terra: Precipitation of Coldest Period
#' @description Calculates precipitation sum of the period with the lowest temperature sum.
#' @param wet spatRaster of precipitation period sums.
#' @param coldest_period Vector indicating the index (1-based) of the coldest period.
#' @return spatRaster with "bio19".
#' @keywords internal
bio19_terra <- function(wet, coldest_period) {
  quihicha_aca <- terra::selectRange(wet, coldest_period)
  names(quihicha_aca) <- "bio19"
  return(quihicha_aca)
}

#' @title bio20_terra: Mean Solar Radiation of Units
#' @description Calculates mean solar radiation across all temporal units.
#' @param srad spatRaster of average solar radiation for each unit.
#' @return spatRaster with "bio20".
#' @keywords internal
bio20_terra <- function(srad) {
  gueta <- terra::app(srad, mean, na.rm = TRUE)
  names(gueta) <- "bio20"
  return(gueta)
}

#' @title bio21_terra: Highest Solar Radiation Unit
#' @description Identifies highest solar radiation unit, potentially using a static index.
#' @param srad spatRaster of solar radiation values for each unit.
#' @return spatRaster with "bio21".
#' @keywords internal
bio21_terra <- function(srad, high_rad_unit = NULL) {
  if (!is.null(high_rad_unit)) {
    gueta_ata <- terra::selectRange(srad, high_rad_unit)
  } else {
    gueta_ata <- terra::app(srad, max, na.rm = TRUE)
  }
  names(gueta_ata) <- "bio21"
  return(gueta_ata)
}

#' @title bio22_terra: Lowest Solar Radiation Unit
#' @description Identifies lowest solar radiation unit, potentially using a static index.
#' @param srad spatRaster of solar radiation values for each unit.
#' @return spatRaster with "bio22".
#' @keywords internal
bio22_terra <- function(srad, low_rad_unit = NULL) {
  if (!is.null(low_rad_unit)) {
    gueta_bosa <- terra::selectRange(srad, low_rad_unit)
  } else {
    gueta_bosa <- terra::app(srad, min, na.rm = TRUE)
  }
  names(gueta_bosa) <- "bio22"
  return(gueta_bosa)
}

#' @title bio23_terra: Solar Radiation Seasonality (CV)
#' @description Calculates coefficient of variation in solar radiation across units.
#' @param srad Matrix containing solar radiation values for each unit.
#' @return spatRaster with "bio23".
#' @keywords internal
bio23_terra <- function(srad) {
  gueta_mica <- cv_cli(srad)
  names(gueta_mica) <- "bio23"
  return(gueta_mica)
}

#' @title bio24_terra: Solar Radiation of Wettest Period
#' @description Calculates solar radiation mean of the period with the highest precipitation sum.
#' @param prad spatRaster of solar radiation period means.
#' @param wettest_period Vector indicating the index (1-based) of the wettest period.
#' @return spatRaster with "bio24".
#' @keywords internal
bio24_terra <- function(prad, wettest_period) {
  gueta_muihica <- terra::selectRange(prad, wettest_period)
  names(gueta_muihica) <- "bio24"
  return(gueta_muihica)
}

#' @title bio25_terra: Solar Radiation of Driest Period
#' @description Calculates solar radiation mean of the period with the highest precipitation sum.
#' @param prad spatRaster of solar radiation period means.
#' @param driest_period Vector indicating the index (1-based) of the driest period.
#' @return spatRaster with "bio25".
#' @keywords internal
bio25_terra <- function(prad, driest_period) {
  gueta_hisca <- terra::selectRange(prad, driest_period)
  names(gueta_hisca) <- "bio25"
  return(gueta_hisca)
}

#' @title bio26_terra: Solar Radiation of Warmest Period
#' @description Calculates solar radiation mean of the period with the highest temperature mean.
#' @param prad spatRaster of solar radiation period means.
#' @param warmest_period Vector indicating the index (1-based) of the warmest period.
#' @return spatRaster with "bio26".
#' @keywords internal
bio26_terra <- function(prad, warmest_period) {
  gueta_ta <- terra::selectRange(prad, warmest_period)
  names(gueta_ta) <- "bio26"
  return(gueta_ta)
}

#' @title bio27_terra: Solar Radiation of Coldest Period
#' @description Calculates solar radiation mean of the period with the lowest temperature mean.
#' @param prad spatRaster of solar radiation period means.
#' @param coldest_period Vector indicating the index (1-based) of the coldest period.
#' @return spatRaster with "bio27".
#' @keywords internal
bio27_terra <- function(prad, coldest_period) {
  gueta_cuhupcua <- terra::selectRange(prad, coldest_period)
  names(gueta_cuhupcua) <- "bio27"
  return(gueta_cuhupcua)
}

#' @title bio28_terra: Mean Moisture of Units
#' @description Calculates mean moisture across all temporal units.
#' @param mois spatRaster of average moisture for each unit.
#' @return spatRaster with "bio28".
#' @keywords internal
bio28_terra <- function(mois) {
  gueta_suhusa <- terra::app(mois, mean, na.rm = TRUE)
  names(gueta_suhusa) <- "bio28"
  return(gueta_suhusa)
}

#' @title bio29_terra: Highest Moisture Unit
#' @description Identifies highest moisture unit, potentially using a static index.
#' @param mois spatRaster of moisture values for each unit.
#' @return spatRaster with "bio29".
#' @keywords internal
bio29_terra <- function(mois, high_mois_unit = NULL) {
  if (!is.null(high_mois_unit)) {
    gueta_aca <- terra::selectRange(mois, high_mois_unit)
  } else {
    gueta_aca <- terra::app(mois, max, na.rm = TRUE)
  }
  names(gueta_aca) <- "bio29"
  return(gueta_aca)
}

#' @title bio30_terra: Lowest Moisture Unit
#' @description Identifies lowest moisture unit, potentially using a static index.
#' @param mois spatRaster of moisture values for each unit.
#' @return spatRaster with "bio30".
#' @keywords internal
bio30_terra <- function(mois, low_mois_unit = NULL) {
  if (!is.null(low_mois_unit)) {
    gueta_ubchihica <- terra::selectRange(mois, low_mois_unit)
  } else {
    gueta_ubchihica <- terra::app(mois, min, na.rm = TRUE)
  }
  names(gueta_ubchihica) <- "bio30"
  return(gueta_ubchihica)
}

#' @title bio31_terra: Moisture Seasonality (Standard Deviation)
#' @description Calculates coefficient of variation in moisture across units.
#' @param mois Matrix containing moisture values for each unit.
#' @return spatRaster with "bio31".
#' @keywords internal
bio31_terra <- function(mois) {
  gueta_quihicha_ata <- 100 * terra::stdev(mois, pop = FALSE, na.rm = TRUE)
  names(gueta_quihicha_ata) <- "bio31"
  return(gueta_quihicha_ata)
}

#' @title bio32_terra: Moisture of the Most Moist Period
#' @description Calculates moisture mean of the most moist period.
#' @param pmois spatRaster of moisture period means.
#' @param high_mois_period Vector indicating the index (1-based) of the most moist period.
#' @return spatRaster with "bio32".
#' @keywords internal
bio32_terra <- function(pmois, high_mois_period) {
  gueta_quihicha_bosa <- terra::selectRange(pmois, high_mois_period)
  names(gueta_quihicha_bosa) <- "bio32"
  return(gueta_quihicha_bosa)
}

#' @title bio33_terra: Moisture of the Least Moist Period
#' @description Calculates moisture mean of the least moist period.
#' @param pmois spatRaster of moisture period means.
#' @param low_mois_period Vector indicating the index (1-based) of the least moist period.
#' @return spatRaster with "bio33".
#' @keywords internal
bio33_terra <- function(pmois, low_mois_period) {
  gueta_quihicha_mica <- terra::selectRange(pmois, low_mois_period)
  names(gueta_quihicha_mica) <- "bio33"
  return(gueta_quihicha_mica)
}

#' @title bio34_terra: Moisture of Warmest Period
#' @description Calculates moisture mean of the period with the highest temperature mean.
#' @param pmois spatRaster of moisture period means.
#' @param warmest_period Vector indicating the index (1-based) of the warmest period.
#' @return spatRaster with "bio34".
#' @keywords internal
bio34_terra <- function(pmois, warmest_period) {
  gueta_quihicha_muihica <- terra::selectRange(pmois, warmest_period)
  names(gueta_quihicha_muihica) <- "bio34"
  return(gueta_quihicha_muihica)
}

#' @title bio35_terra: Moisture of Coldest Period
#' @description Calculates moisture mean of the period with the lowest temperature mean.
#' @param pmois spatRaster of moisture period means.
#' @param coldest_period Vector indicating the index (1-based) of the coldest period.
#' @return spatRaster with "bio35".
#' @keywords internal
bio35_terra <- function(pmois, coldest_period) {
  gueta_quihicha_hisca <- terra::selectRange(pmois, coldest_period)
  names(gueta_quihicha_hisca) <- "bio35"
  return(gueta_quihicha_hisca)
}