% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bios_fast.R
\name{bio04_fast}
\alias{bio04_fast}
\title{bio04_fast: Temperature Seasonality (Std Dev * 100)}
\usage{
bio04_fast(tavg, cell)
}
\arguments{
\item{tavg}{Matrix of average temperatures for each unit.}

\item{cell}{Vector of original cell IDs.}
}
\value{
Matrix with "bio04", "cell".
}
\description{
Calculates the standard deviation of average temperatures across units, multiplied by 100.
}
\keyword{internal}
