% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bios_fast.R
\name{bio20_fast}
\alias{bio20_fast}
\title{bio20_fast: Mean Solar Radiation of Units}
\usage{
bio20_fast(srad, cell)
}
\arguments{
\item{srad}{Matrix of average solar radiation for each unit.}

\item{cell}{Vector of original cell IDs.}
}
\value{
Matrix with columns: "bio20", "cell".
}
\description{
Calculates mean solar radiation across all temporal units.
}
\keyword{internal}
