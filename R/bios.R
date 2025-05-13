# Function to create analogous bioclimatic variables
# ----------
# P1. Mean Temperature of Units
#' @export
bio_01 <- function(tavg) {
  ata <- terra::app(tavg, mean, na.rm = TRUE)
  names(ata) <- "bio01"
  return(ata)
}

# P2. Mean Diurnal Range(Mean(unit max-min))
#' @export
bio_02 <- function(tmin, tmax) {
  bosa <- terra::app(tmax - tmin, mean, na.rm = TRUE)
  names(bosa) <- "bio02"
  return(bosa)
}

# P3. Isothermality (P2 / P7)
#' @export
bio_03 <- function(bio02, bio07) {
   mica <- 100 * bio02 / bio07
  names(mica) <- "bio03"
  return(mica)
}

# P4. Temperature Seasonality (standard deviation)
#' @export
bio_04 <- function(tavg) {
  muihica <- 100 * terra::stdev(tavg, pop = FALSE, na.rm = TRUE)
  names(muihica) <- "bio04"
  return(muihica)
}

# P5. Max Temperature of Warmest Unit
#' @export
bio_05 <- function(tmax, warmest_unit = NULL) {
  if (!is.null(warmest_unit)) {
    hisca <- terra::selectRange(tmax, warmest_unit)
  } else {
    hisca <- terra::app(tmax, max, na.rm = TRUE)
  }
  names(hisca) <- "bio05"
  return(hisca)
}

# P6. Min Temperature of Coldest Unit
#' @export
bio_06 <- function(tmin, coldest_unit = NULL) {
  if (!is.null(coldest_unit)) {
    ta <- terra::selectRange(tmin, coldest_unit)
  } else {
    ta <- terra::app(tmin, min, na.rm = TRUE)
  }
  names(ta) <- "bio06"
  return(ta)
}

# P7. Temperature Range of Units (P5-P6)
#' @export
bio_07 <- function(bio05, bio06) {
  cuhupcua <- bio05 - bio06
  names(cuhupcua) <- "bio07"
  return(cuhupcua)
}

# P08. Mean Temperature of Wettest Period
#' @export
bio_08 <- function(tmp, wettest_period) {
  suhusa <- terra::selectRange(tmp, wettest_period)
  names(suhusa) <- "bio08"
  return(suhusa)
}

# P09. Mean Temperature of Driest Period
#' @export
bio_09 <- function(tmp, driest_period) {
  aca <- terra::selectRange(tmp, driest_period)
  names(aca) <- "bio09"
  return(aca)
}

# P10. Mean Temperature of Warmest Period
#' @export
bio_10 <- function(tmp, warmest_period) {
  ubchihica <- terra::selectRange(tmp, warmest_period)
  names(ubchihica) <- "bio10"
  return(ubchihica)
}

# P11. Mean Temperature of Coldest Period
#' @export
bio_11 <- function(tmp, coldest_period) {
  quihicha_ata <- terra::selectRange(tmp, coldest_period)
  names(quihicha_ata) <- "bio11"
  return(quihicha_ata)
}

# P12. Precipitation Sum
#' @export
bio_12 <- function(prcp) {
  quihicha_bosa <- terra::app(prcp, sum, na.rm = TRUE)
  names(quihicha_bosa) <- "bio12"
  return(quihicha_bosa)
}

# P13. Precipitation of Wettest Unit
#' @export
bio_13 <- function(prcp, wettest_unit = NULL) {
  if (!is.null(wettest_unit)) {
    quihicha_mica <- terra::selectRange(prcp, wettest_unit)
  } else {
    quihicha_mica <- terra::app(prcp, max, na.rm = TRUE)
  }
  names(quihicha_mica) <- "bio13"
  return(quihicha_mica)
}

# P14. Precipitation of Driest Unit
#' @export
bio_14 <- function(prcp, driest_unit = NULL) {
  if (!is.null(driest_unit)) {
    quihicha_muihica <- terra::selectRange(prcp, driest_unit)
  } else {
    quihicha_muihica <- terra::app(prcp, min, na.rm = TRUE)
  }
  names(quihicha_muihica) <- "bio14"
  return(quihicha_muihica)
}

# P15. Precipitation Seasonality (Coefficient of Variation)
# the "1 +" is to avoid strange CVs for areas where mean rainfaill is < 1)
#' @export
bio_15 <- function(prcp) {
  quihicha_hisca <- cv_cli(prcp)
  names(quihicha_hisca) <- "bio15"
  return(quihicha_hisca)
}

# P16. Precipitation of Wettest Period
#' @export
bio_16 <- function(wet, wettest_period) {
  quihicha_ta <- terra::selectRange(wet, wettest_period)
  names(quihicha_ta) <- "bio16"
  return(quihicha_ta)
}

# P17. Precipitation of Driest Period
#' @export
bio_17 <- function(wet, driest_period) {
  quihicha_cuhupcua <- terra::selectRange(wet, driest_period)
  names(quihicha_cuhupcua) <- "bio17"
  return(quihicha_cuhupcua)
}

# P18. Precipitation of Warmest Period
#' @export
bio_18 <- function(wet, warmest_period) {
  quihicha_suhusa <- terra::selectRange(wet, warmest_period)
  names(quihicha_suhusa) <- "bio18"
  return(quihicha_suhusa)
}

# P19. Precipitation of Coldest Period
#' @export
bio_19 <- function(wet, coldest_period) {
  quihicha_aca <- terra::selectRange(wet, coldest_period)
  names(quihicha_aca) <- "bio19"
  return(quihicha_aca)
}

# P20. Mean Radiation of Units
#' @export
bio_20 <- function(srad) {
  gueta <- terra::app(srad, mean, na.rm = TRUE)
  names(gueta) <- "bio20"
  return(gueta)
}

# P21. Highest Radiation Unit
#' @export
bio_21 <- function(srad, high_rad_unit = NULL) {
  if (!is.null(high_rad_unit)) {
    gueta_ata <- terra::selectRange(srad, high_rad_unit)
  } else {
    gueta_ata <- terra::app(srad, max, na.rm = TRUE)
  }
  names(gueta_ata) <- "bio21"
  return(gueta_ata)
}

# P22. Lowest Radiation Unit
#' @export
bio_22 <- function(srad, low_rad_unit = NULL) {
  if (!is.null(low_rad_unit)) {
    gueta_bosa <- terra::selectRange(srad, low_rad_unit)
  } else {
    gueta_bosa <- terra::app(srad, min, na.rm = TRUE)
  }
  names(gueta_bosa) <- "bio22"
  return(gueta_bosa)
}

# P23. Radiation Seasonality (Coefficient of Variation)
#' @export
bio_23 <- function(srad) {
  gueta_mica <- cv_cli(srad)
  names(gueta_mica) <- "bio23"
  return(gueta_mica)
}

# P24. Radiation of Wettest Period
#' @export
bio_24 <- function(prad, wettest_period) {
  gueta_muihica <- terra::selectRange(prad, wettest_period)
  names(gueta_muihica) <- "bio24"
  return(gueta_muihica)
}

# P25. Radiation of Driest Period
#' @export
bio_25 <- function(prad, driest_period) {
  gueta_hisca <- terra::selectRange(prad, driest_period)
  names(gueta_hisca) <- "bio25"
  return(gueta_hisca)
}

# P26. Radiation of Warmest Period
#' @export
bio_26 <- function(prad, warmest_period) {
  gueta_ta <- terra::selectRange(prad, warmest_period)
  names(gueta_ta) <- "bio26"
  return(gueta_ta)
}

# P27. Radiation of Coldest Period
#' @export
bio_27 <- function(prad, coldest_period) {
  gueta_cuhupcua <- terra::selectRange(prad, coldest_period)
  names(gueta_cuhupcua) <- "bio27"
  return(gueta_cuhupcua)
}

# P28. Mean Moisture
#' @export
bio_28 <- function(mois) {
  gueta_suhusa <- terra::app(mois, mean, na.rm = TRUE)
  names(gueta_suhusa) <- "bio28"
  return(gueta_suhusa)
}

# P29. Highest Moisture Unit
#' @export
bio_29 <- function(mois, high_mois_unit = NULL) {
  if (!is.null(high_mois_unit)) {
    gueta_aca <- terra::selectRange(mois, high_mois_unit)
  } else {
    gueta_aca <- terra::app(mois, max, na.rm = TRUE)
  }
  names(gueta_aca) <- "bio29"
  return(gueta_aca)
}

# P30. Lowest Moisture Unit
#' @export
bio_30 <- function(mois, low_mois_unit = NULL) {
  if (!is.null(low_mois_unit)) {
    gueta_ubchihica <- terra::selectRange(mois, low_mois_unit)
  } else {
    gueta_ubchihica <- terra::app(mois, min, na.rm = TRUE)
  }
  names(gueta_ubchihica) <- "bio30"
  return(gueta_ubchihica)
}

# P31. Moisture Seasonality (Coefficient of Variation)
#' @export
bio_31 <- function(mois) {
  gueta_quihicha_ata <- cv_cli(mois)
  names(gueta_quihicha_ata) <- "bio31"
  return(gueta_quihicha_ata)
}

# P32. Mean Moisture of the Most Moist Period
#' @export
bio_32 <- function(pmois, high_mois_period) {
  gueta_quihicha_bosa <- terra::selectRange(pmois, high_mois_period)
  names(gueta_quihicha_bosa) <- "bio32"
  return(gueta_quihicha_bosa)
}

# P33. Mean Moisture of the Least Moist Period
#' @export
bio_33 <- function(pmois, low_mois_period) {
  gueta_quihicha_mica <- terra::selectRange(pmois, low_mois_period)
  names(gueta_quihicha_mica) <- "bio33"
  return(gueta_quihicha_mica)
}

# P34. Mean Moisture of Warmest Period
#' @export
bio_34 <- function(pmois, warmest_period) {
  gueta_quihicha_muihica <- terra::selectRange(pmois, warmest_period)
  names(gueta_quihicha_muihica) <- "bio34"
  return(gueta_quihicha_muihica)
}

# P35. Mean Moisture of Coldest Period
#' @export
bio_35 <- function(pmois, coldest_period) {
  gueta_quihicha_hisca <- terra::selectRange(pmois, coldest_period)
  names(gueta_quihicha_hisca) <- "bio35"
  return(gueta_quihicha_hisca)
}