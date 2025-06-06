% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help2.R
\name{bio02_fun}
\alias{bio02_fun}
\title{bio02: Mean Diurnal Range}
\usage{
bio02_fun(tmin, tmax, cell)
}
\arguments{
\item{tmin}{Matrix of minimum temperatures for each unit.}

\item{tmax}{Matrix of maximum temperatures for each unit.}

\item{cell}{Vector of original cell IDs.}
}
\value{
Matrix with "bio02", "cell".
}
\description{
Calculates the mean of (tmax - tmin) across all temporal units.
}
\keyword{internal}
