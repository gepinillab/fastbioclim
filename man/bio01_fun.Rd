% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help2.R
\name{bio01_fun}
\alias{bio01_fun}
\title{bio01: Mean Temperature of Units}
\usage{
bio01_fun(tavg, cell)
}
\arguments{
\item{tavg}{Matrix of average temperatures for each unit.}

\item{cell}{Vector of original cell IDs.}
}
\value{
Matrix with columns: "bio01", "cell".
}
\description{
Calculates mean temperature across all temporal units.
}
\keyword{internal}
