% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_terra.R
\name{stats_terra}
\alias{stats_terra}
\title{In-Memory Custom Variable Summarization}
\usage{
stats_terra(
  variable,
  stats = c("mean", "max", "min", "cv_cli", "max_period", "min_period"),
  period_length = 3,
  period_stats = "mean",
  circular = TRUE,
  inter_variable = NULL,
  inter_stats = c("max_inter", "min_inter"),
  max_unit = NULL,
  min_unit = NULL,
  max_period = NULL,
  min_period = NULL,
  max_interactive = NULL,
  min_interactive = NULL,
  prefix_variable = "var",
  suffix_inter_max = "inter_high",
  suffix_inter_min = "inter_low",
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE,
  output_dir = tempdir(),
  ...
)
}
\arguments{
\item{variable}{A SpatRaster for the primary variable.}

\item{stats}{Character vector of stats to compute.}

\item{prefix_variable}{Character, prefix for output layer names.}

\item{...}{Other arguments including inter_variable, period_length, circular,
static index SpatRasters, etc.}
}
\value{
A `terra::SpatRaster` object containing the calculated summary layers.
}
\description{
Internal function to calculate custom summary statistics using `terra` functions.
Designed for datasets that fit into RAM.
}
\seealso{
The user-facing wrapper function `derived_statistics()`.
}
\keyword{internal}
