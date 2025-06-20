% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bios_fast.R
\name{bio24_fast}
\alias{bio24_fast}
\title{bio24_fast: Solar Radiation of Wettest Period}
\usage{
bio24_fast(speriod, pperiod_max_idx, cell)
}
\arguments{
\item{speriod}{Matrix of solar radiation period means (output from `var_periods`).}

\item{pperiod_max_idx}{Vector indicating the index (1-based) of the wettest period.}

\item{cell}{Vector of original cell IDs.}
}
\value{
Matrix with "bio24", "cell".
}
\description{
Calculates solar radiation mean of the period with the highest precipitation sum.
}
\keyword{internal}
