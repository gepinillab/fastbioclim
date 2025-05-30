% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help2.R
\name{bio19_fun}
\alias{bio19_fun}
\title{bio19: Precipitation of Coldest Period}
\usage{
bio19_fun(pperiod, tperiod_min_idx, cell)
}
\arguments{
\item{pperiod}{Matrix of precipitation period sums.}

\item{tperiod_min_idx}{Vector indicating the index (1-based) of the coldest period.}

\item{cell}{Vector of original cell IDs.}
}
\value{
Matrix with "bio19", "cell".
}
\description{
Calculates precipitation sum of the period with the lowest temperature sum.
}
\keyword{internal}
