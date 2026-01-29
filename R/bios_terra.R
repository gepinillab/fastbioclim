# Function to create analogous bioclimatic variables
# ----------
#' @title bio01_terra: Mean Temperature
#' @description Calculates mean temperature across all temporal units (layers).
#' @param tavg A `SpatRaster` object where each layer represents average temperature for a temporal unit (e.g., 12 months).
#' @return A single-layer `SpatRaster` with the calculated mean temperature, named "bio01".
#' @keywords internal
bio01_terra <- function(tavg) {
  ata <- terra::app(tavg, mean, na.rm = TRUE)
  names(ata) <- "bio01"
  return(ata)
}

#' @title bio02_terra: Mean Diurnal Range
#' @description Calculates the mean of the diurnal temperature range (tmax - tmin) across all temporal units (layers).
#' @param tmin A `SpatRaster` object of minimum temperatures, where each layer is a temporal unit. 
#'    Must have the same dimensions and number of layers as `tmax`.
#' @param tmax A `SpatRaster` object of maximum temperatures, where each layer is a temporal unit. 
#'    Must have the same dimensions and number of layers as `tmin`.
#' @return A single-layer `SpatRaster` with the mean diurnal range, named "bio02".
#' @keywords internal
bio02_terra <- function(tmin, tmax) {
  bosa <- terra::app(tmax - tmin, mean, na.rm = TRUE)
  names(bosa) <- "bio02"
  return(bosa)
}

#' @title bio03_terra: Isothermality
#' @description Calculates Isothermality, defined as (bio02 / bio07) * 100.
#' @param bio02 A single-layer `SpatRaster` of Mean Diurnal Range (Bio02).
#' @param bio07 A single-layer `SpatRaster` of Temperature Range (Bio07).
#' @return A single-layer `SpatRaster` with the calculated isothermality, named "bio03".
#' @keywords internal
bio03_terra <- function(bio02, bio07) {
   mica <- 100 * bio02 / bio07
  names(mica) <- "bio03"
  return(mica)
}

#' @title bio04_terra: Temperature Seasonality (Std Dev * 100)
#' @description Calculates the standard deviation of average temperatures across all layers, multiplied by 100.
#' @param tavg A `SpatRaster` object where each layer represents average temperature for a temporal unit.
#' @return A single-layer `SpatRaster` with the temperature seasonality, named "bio04".
#' @keywords internal
bio04_terra <- function(tavg) {
  muihica <- 100 * terra::stdev(tavg, pop = FALSE, na.rm = TRUE)
  names(muihica) <- "bio04"
  return(muihica)
}

#' @title bio05_terra: Max Temperature of Warmest Unit
#' @description Identifies the maximum temperature of the warmest temporal unit (layer).
#' @details This function calculates bio05 following the standard definition used by WorldClim 
#'   (Hijmans et al., 2005) and ANUCLIM 6.1 (Xu & Hutchinson, 2013), which is the single highest value from all 
#'   maximum temperature layers (a "max of maxes"). It does not use mean 
#'   temperature to first identify the warmest month.
#' @param tmax A `SpatRaster` object of maximum temperatures, where each layer is a temporal unit.
#' @param warmest_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a 
#'    static layer index (1-based) from which to extract the value. If `NULL` (the default), the overall maximum across all 
#'    layers is calculated.
#' @return A single-layer `SpatRaster` with the maximum temperature, named "bio05".
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
#' @description Identifies the minimum temperature of the coldest temporal unit (layer).
#' @details This function calculates bio06 following the standard definition used by WorldClim 
#'   (Hijmans et al., 2005; O'Donnell & Ignizio, 2012) and ANUCLIM 6.1 (Xu & Hutchinson, 2013), which is the 
#'   single lowest value from all 
#'   minimum temperature layers (a "min of mins"). It does not use mean 
#'   temperature to first identify the coldest month.
#' @param tmin A `SpatRaster` object of minimum temperatures, where each layer is a temporal unit.
#' @param coldest_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a static layer index (1-based) from which to extract the value. If `NULL` (the default), the overall minimum across all layers is calculated.
#' @return A single-layer `SpatRaster` with the minimum temperature, named "bio06".
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
bio06_terra <- function(tmin, coldest_unit = NULL) {
  if (!is.null(coldest_unit)) {
    ta <- terra::selectRange(tmin, coldest_unit)
  } else {
    ta <- terra::app(tmin, min, na.rm = TRUE)
  }
  names(ta) <- "bio06"
  return(ta)
}

#' @title bio07_terra: Temperature Range (bio05 - bio06)
#' @description Calculates the difference between the Maximum Temperature of the Warmest Unit (bio05) and the Minimum Temperature of the Coldest Unit (bio06).
#' @param bio05 A single-layer `SpatRaster` of Bio05 values.
#' @param bio06 A single-layer `SpatRaster` of Bio06 values.
#' @return A single-layer `SpatRaster` with the temperature annual range, named "bio07".
#' @keywords internal
bio07_terra <- function(bio05, bio06) {
  cuhupcua <- bio05 - bio06
  names(cuhupcua) <- "bio07"
  return(cuhupcua)
}

#' @title bio08_terra: Mean Temperature of Wettest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the wettest.
#' @param tmp A `SpatRaster` object where each layer represents the mean temperature for a rolling period.
#' @param wettest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the wettest period.
#' @return A single-layer `SpatRaster` containing the mean temperature of the wettest period, named "bio08".
#' @keywords internal
bio08_terra <- function(tmp, wettest_period) {
  suhusa <- terra::selectRange(tmp, wettest_period)
  names(suhusa) <- "bio08"
  return(suhusa)
}

#' @title bio09_terra: Mean Temperature of Driest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the driest.
#' @param tmp A `SpatRaster` object where each layer represents the mean temperature for a rolling period.
#' @param driest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the driest period.
#' @return A single-layer `SpatRaster` containing the mean temperature of the driest period, named "bio09".
#' @keywords internal
bio09_terra <- function(tmp, driest_period) {
  aca <- terra::selectRange(tmp, driest_period)
  names(aca) <- "bio09"
  return(aca)
}

#' @title bio10_terra: Mean Temperature of Warmest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the warmest.
#' @param tmp A `SpatRaster` object where each layer represents the mean temperature for a rolling period.
#' @param warmest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the warmest period.
#' @return A single-layer `SpatRaster` containing the mean temperature of the warmest period, named "bio10".
#' @keywords internal
bio10_terra <- function(tmp, warmest_period) {
  ubchihica <- terra::selectRange(tmp, warmest_period)
  names(ubchihica) <- "bio10"
  return(ubchihica)
}

#' @title bio11_terra: Mean Temperature of Coldest Period
#' @description Calculates the mean temperature of the specific rolling period identified as the coldest.
#' @param tmp A `SpatRaster` object where each layer represents the mean temperature for a rolling period.
#' @param coldest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the coldest period.
#' @return A single-layer `SpatRaster` containing the mean temperature of the coldest period, named "bio11".
#' @keywords internal
bio11_terra <- function(tmp, coldest_period) {
  quihicha_ata <- terra::selectRange(tmp, coldest_period)
  names(quihicha_ata) <- "bio11"
  return(quihicha_ata)
}

#' @title bio12_terra: Total Precipitation
#' @description Calculates the total precipitation (sum) across all temporal units (layers).
#' @param prcp A `SpatRaster` object where each layer represents precipitation for a temporal unit.
#' @return A single-layer `SpatRaster` with the total precipitation, named "bio12".
#' @keywords internal
bio12_terra <- function(prcp) {
  quihicha_bosa <- terra::app(prcp, sum, na.rm = TRUE)
  names(quihicha_bosa) <- "bio12"
  return(quihicha_bosa)
}

#' @title bio13_terra: Precipitation of Wettest Unit
#' @description Identifies the precipitation of the wettest temporal unit (layer).
#' @param prcp A `SpatRaster` object where each layer represents precipitation for a temporal unit.
#' @param wettest_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a static layer index (1-based) from which to extract the value. If `NULL` (the default), the overall maximum across all layers is calculated.
#' @return A single-layer `SpatRaster` with the precipitation of the wettest unit, named "bio13".
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
#' @description Identifies the precipitation of the driest temporal unit (layer).
#' @param prcp A `SpatRaster` object where each layer represents precipitation for a temporal unit.
#' @param driest_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a static layer index (1-based) from which to extract the value. If `NULL` (the default), the overall minimum across all layers is calculated.
#' @return A single-layer `SpatRaster` with the precipitation of the driest unit, named "bio14".
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
#' @description Calculates the Coefficient of Variation (CV) of precipitation across all layers.
#' @param prcp A `SpatRaster` object where each layer represents precipitation for a temporal unit.
#' @note The formula adds 1 to the mean to avoid division by zero in arid areas.
#' @return A single-layer `SpatRaster` with the precipitation seasonality, named "bio15".
#' @keywords internal
bio15_terra <- function(prcp) {
  quihicha_hisca <- cv_cli(prcp)
  names(quihicha_hisca) <- "bio15"
  return(quihicha_hisca)
}

#' @title bio16_terra: Precipitation of Wettest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the wettest.
#' @param wet A `SpatRaster` object where each layer is the precipitation sum for a rolling period.
#' @param wettest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the wettest period.
#' @return A single-layer `SpatRaster` with the precipitation of the wettest period, named "bio16".
#' @keywords internal
bio16_terra <- function(wet, wettest_period) {
  quihicha_ta <- terra::selectRange(wet, wettest_period)
  names(quihicha_ta) <- "bio16"
  return(quihicha_ta)
}

#' @title bio17_terra: Precipitation of Driest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the driest.
#' @param wet A `SpatRaster` object where each layer is the precipitation sum for a rolling period.
#' @param driest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the driest period.
#' @return A single-layer `SpatRaster` with the precipitation of the driest period, named "bio17".
#' @keywords internal
bio17_terra <- function(wet, driest_period) {
  quihicha_cuhupcua <- terra::selectRange(wet, driest_period)
  names(quihicha_cuhupcua) <- "bio17"
  return(quihicha_cuhupcua)
}

#' @title bio18_terra: Precipitation of Warmest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the warmest.
#' @param wet A `SpatRaster` object where each layer is the precipitation sum for a rolling period.
#' @param warmest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the warmest period.
#' @return A single-layer `SpatRaster` with the precipitation of the warmest period, named "bio18".
#' @keywords internal
bio18_terra <- function(wet, warmest_period) {
  quihicha_suhusa <- terra::selectRange(wet, warmest_period)
  names(quihicha_suhusa) <- "bio18"
  return(quihicha_suhusa)
}

#' @title bio19_terra: Precipitation of Coldest Period
#' @description Calculates the total precipitation of the specific rolling period identified as the coldest.
#' @param wet A `SpatRaster` object where each layer is the precipitation sum for a rolling period.
#' @param coldest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the coldest period.
#' @return A single-layer `SpatRaster` with the precipitation of the coldest period, named "bio19".
#' @keywords internal
bio19_terra <- function(wet, coldest_period) {
  quihicha_aca <- terra::selectRange(wet, coldest_period)
  names(quihicha_aca) <- "bio19"
  return(quihicha_aca)
}

#' @title bio20_terra: Mean Radiation
#' @description Calculates mean solar radiation across all temporal units (layers).
#' @param srad A `SpatRaster` object where each layer represents solar radiation for a temporal unit.
#' @return A single-layer `SpatRaster` with the mean solar radiation, named "bio20".
#' @keywords internal
bio20_terra <- function(srad) {
  gueta <- terra::app(srad, mean, na.rm = TRUE)
  names(gueta) <- "bio20"
  return(gueta)
}

#' @title bio21_terra: Highest Moisture Unit
#' @description Identifies the highest solar radiation of the unit with the highest value.
#' @param srad A `SpatRaster` object where each layer represents solar radiation for a temporal unit.
#' @param high_rad_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a static layer index (1-based). If `NULL`, the overall maximum across all layers is calculated.
#' @return A single-layer `SpatRaster` with the maximum solar radiation, named "bio21".
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

#' @title bio22_terra: Lowest Radiation Unit
#' @description Identifies the lowest solar radiation of the unit with the lowest value.
#' @param srad A `SpatRaster` object where each layer represents solar radiation for a temporal unit.
#' @param low_rad_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a static layer index (1-based). If `NULL`, the overall minimum across all layers is calculated.
#' @return A single-layer `SpatRaster` with the minimum solar radiation, named "bio22".
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

#' @title bio23_terra: Radiation Seasonality (CV)
#' @description Calculates the Coefficient of Variation (CV) of solar radiation across all layers.
#' @param srad A `SpatRaster` object where each layer represents solar radiation for a temporal unit.
#' @return A single-layer `SpatRaster` with the solar radiation seasonality, named "bio23".
#' @keywords internal
bio23_terra <- function(srad) {
  gueta_mica <- cv_cli(srad)
  names(gueta_mica) <- "bio23"
  return(gueta_mica)
}

#' @title bio24_terra: Radiation of Wettest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the wettest.
#' @param prad A `SpatRaster` object where each layer is the mean solar radiation for a rolling period.
#' @param wettest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the wettest period.
#' @return A single-layer `SpatRaster` with the solar radiation of the wettest period, named "bio24".
#' @keywords internal
bio24_terra <- function(prad, wettest_period) {
  gueta_muihica <- terra::selectRange(prad, wettest_period)
  names(gueta_muihica) <- "bio24"
  return(gueta_muihica)
}

#' @title bio25_terra: Radiation of Driest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the driest.
#' @param prad A `SpatRaster` object where each layer is the mean solar radiation for a rolling period.
#' @param driest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the driest period.
#' @return A single-layer `SpatRaster` with the solar radiation of the driest period, named "bio25".
#' @keywords internal
bio25_terra <- function(prad, driest_period) {
  gueta_hisca <- terra::selectRange(prad, driest_period)
  names(gueta_hisca) <- "bio25"
  return(gueta_hisca)
}

#' @title bio26_terra: Radiation of Warmest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the warmest.
#' @param prad A `SpatRaster` object where each layer is the mean solar radiation for a rolling period.
#' @param warmest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the warmest period.
#' @return A single-layer `SpatRaster` with the solar radiation of the warmest period, named "bio26".
#' @keywords internal
bio26_terra <- function(prad, warmest_period) {
  gueta_ta <- terra::selectRange(prad, warmest_period)
  names(gueta_ta) <- "bio26"
  return(gueta_ta)
}

#' @title bio27_terra: Radiation of Coldest Period
#' @description Calculates the mean solar radiation of the specific rolling period identified as the coldest.
#' @param prad A `SpatRaster` object where each layer is the mean solar radiation for a rolling period.
#' @param coldest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the coldest period.
#' @return A single-layer `SpatRaster` with the solar radiation of the coldest period, named "bio27".
#' @keywords internal
bio27_terra <- function(prad, coldest_period) {
  gueta_cuhupcua <- terra::selectRange(prad, coldest_period)
  names(gueta_cuhupcua) <- "bio27"
  return(gueta_cuhupcua)
}

#' @title bio28_terra: Mean Moisture
#' @description Calculates mean moisture across all temporal units (layers).
#' @param mois A `SpatRaster` object where each layer represents moisture for a temporal unit.
#' @return A single-layer `SpatRaster` with the mean moisture, named "bio28".
#' @keywords internal
bio28_terra <- function(mois) {
  gueta_suhusa <- terra::app(mois, mean, na.rm = TRUE)
  names(gueta_suhusa) <- "bio28"
  return(gueta_suhusa)
}

#' @title bio29_terra: Highest Moisture Unit
#' @description Identifies the highest moisture of the unit with the highest value.
#' @param mois A `SpatRaster` object where each layer represents moisture for a temporal unit.
#' @param high_mois_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a static layer index (1-based). If `NULL`, the overall maximum across all layers is calculated.
#' @return A single-layer `SpatRaster` with the maximum moisture, named "bio29".
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
#' @description Identifies the lowest moisture of the unit with the lowest value.
#' @param mois A `SpatRaster` object where each layer represents moisture for a temporal unit.
#' @param low_mois_unit (Optional) A single-layer `SpatRaster` where cell values are integers indicating a static layer index (1-based). If `NULL`, the overall minimum across all layers is calculated.
#' @return A single-layer `SpatRaster` with the minimum moisture, named "bio30".
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

#' @title bio31_terra: Moisture Seasonality (Std Dev * 100)
#' @description Calculates the standard deviation of moisture across all layers, multiplied by 100.
#' @param mois A `SpatRaster` object where each layer represents moisture for a temporal unit.
#' @return A single-layer `SpatRaster` with the moisture seasonality, named "bio31".
#' @keywords internal
bio31_terra <- function(mois) {
  gueta_quihicha_ata <- 100 * terra::stdev(mois, pop = FALSE, na.rm = TRUE)
  names(gueta_quihicha_ata) <- "bio31"
  return(gueta_quihicha_ata)
}

#' @title bio32_terra: Mean Moisture of Most Moist Period
#' @description Calculates the mean moisture of the specific rolling period identified as the most moist.
#' @param pmois A `SpatRaster` object where each layer is the mean moisture for a rolling period.
#' @param high_mois_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the most moist period.
#' @return A single-layer `SpatRaster` with the moisture of the most moist period, named "bio32".
#' @keywords internal
bio32_terra <- function(pmois, high_mois_period) {
  gueta_quihicha_bosa <- terra::selectRange(pmois, high_mois_period)
  names(gueta_quihicha_bosa) <- "bio32"
  return(gueta_quihicha_bosa)
}

#' @title bio33_terra: Mean Moisture of Least Moist Period
#' @description Calculates the mean moisture of the specific rolling period identified as the least moist.
#' @param pmois A `SpatRaster` object where each layer is the mean moisture for a rolling period.
#' @param low_mois_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the least moist period.
#' @return A single-layer `SpatRaster` with the moisture of the least moist period, named "bio33".
#' @keywords internal
bio33_terra <- function(pmois, low_mois_period) {
  gueta_quihicha_mica <- terra::selectRange(pmois, low_mois_period)
  names(gueta_quihicha_mica) <- "bio33"
  return(gueta_quihicha_mica)
}

#' @title bio34_terra: Mean Moisture of Warmest Period
#' @description Calculates the mean moisture of the specific rolling period identified as the warmest.
#' @param pmois A `SpatRaster` object where each layer is the mean moisture for a rolling period.
#' @param warmest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the warmest period.
#' @return A single-layer `SpatRaster` with the moisture of the warmest period, named "bio34".
#' @keywords internal
bio34_terra <- function(pmois, warmest_period) {
  gueta_quihicha_muihica <- terra::selectRange(pmois, warmest_period)
  names(gueta_quihicha_muihica) <- "bio34"
  return(gueta_quihicha_muihica)
}

#' @title bio35_terra: Mean Moisture of Coldest Period
#' @description Calculates the mean moisture of the specific rolling period identified as the coldest.
#' @param pmois A `SpatRaster` object where each layer is the mean moisture for a rolling period.
#' @param coldest_period A single-layer `SpatRaster` where cell values are integers indicating the layer index (1-based) of the coldest period.
#' @return A single-layer `SpatRaster` with the moisture of the coldest period, named "bio35".
#' @keywords internal
bio35_terra <- function(pmois, coldest_period) {
  gueta_quihicha_hisca <- terra::selectRange(pmois, coldest_period)
  names(gueta_quihicha_hisca) <- "bio35"
  return(gueta_quihicha_hisca)
}