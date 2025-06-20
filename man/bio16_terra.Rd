% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bios_terra.R
\name{bio16_terra}
\alias{bio16_terra}
\title{bio16_terra: Precipitation of Wettest Period}
\usage{
bio16_terra(wet, wettest_period)
}
\arguments{
\item{wet}{spatRaster of precipitation period sums.}

\item{wettest_period}{Vector indicating the index (1-based) of the wettest period.}
}
\value{
spatRaster with "bio16".
}
\description{
Calculates precipitation sum of the period with the highest precipitation sum.
}
\keyword{internal}
