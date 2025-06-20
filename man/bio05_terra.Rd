% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bios_terra.R
\name{bio05_terra}
\alias{bio05_terra}
\title{bio05_terra: Max Temperature of Warmest Unit}
\usage{
bio05_terra(tmax, warmest_unit = NULL)
}
\arguments{
\item{tmax}{spatRaster of maximum temperatures for each unit.}
}
\value{
spatRaster with "bio05".
}
\description{
Identifies max temperature of the warmest unit, potentially using a static index.
}
\keyword{internal}
