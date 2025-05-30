% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioclim.R
\name{bioclim_vars}
\alias{bioclim_vars}
\title{Compute Specified Bioclimatic Variables}
\usage{
bioclim_vars(
  bios,
  n_units,
  tmin_path = NULL,
  tmax_path = NULL,
  prec_path = NULL,
  tavg_path = NULL,
  srad_path = NULL,
  mois_path = NULL,
  period_length = 3,
  circular = TRUE,
  user_region = NULL,
  tile_degrees = 5,
  temp_dir = tempdir(),
  write_raw_vars = FALSE,
  ...
)
}
\arguments{
\item{bios}{Numeric vector specifying which bioclimatic variables (1-35) to compute.}

\item{n_units}{Integer. The number of temporal units (layers) per input variable.}

\item{tmin_path}{Character vector of `n_units` paths to minimum temperature rasters.}

\item{tmax_path}{Character vector of `n_units` paths to maximum temperature rasters.}

\item{prec_path}{Character vector of `n_units` paths to precipitation rasters.}

\item{tavg_path}{Character vector of `n_units` paths to optional average temperature rasters.}

\item{srad_path}{Character vector of `n_units` paths to solar radiation rasters.}

\item{mois_path}{Character vector of `n_units` paths to moisture rasters.}

\item{period_length}{Integer. The number of units defining a "period". Default: 3.}

\item{circular}{Logical. Calculate periods wrapping around the cycle? Default: TRUE.}

\item{user_region}{Optional. An `sf` or `terra::SpatVector` object defining the
processing area. If `NULL`, the full extent of the input rasters is used.}

\item{tile_degrees}{Numeric. Approximate size of processing tiles in degrees. Default: 5.}

\item{temp_dir}{Character. Path for temporary tile files. Default: `tempdir()`.}

\item{write_raw_vars}{Logical. Save intermediate extracted climate data. Default: FALSE.}

\item{...}{Additional arguments, including optional paths to static index rasters
(e.g., `warmest_period_path`).}
}
\value{
Character string: Path to the temporary directory containing intermediate `.qs` files.
}
\description{
Calculates specified bioclimatic variables (1-35) from climate rasters
with user-defined temporal units, optionally using static indices. AOI is defined
by `user_region` or defaults to the full raster extent. Uses parallel processing,
spatial tiling, and `exactextractr`.
}
\details{
Calculates BIOs 1-35. The Area of Interest (AOI) is determined by `user_region`.
  If `user_region` is NULL, the AOI is the full spatial extent of the input climate rasters.
  Static indices provided via `...` override dynamic calculations.
  Uses `exactextractr` for extraction and `Rfast` for matrix calculations.
  Assemble results using a compatible `write_layers` function.
}
