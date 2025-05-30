% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help2.R
\name{bio03_fun}
\alias{bio03_fun}
\title{bio03: Isothermality}
\usage{
bio03_fun(bio02V, bio07V, cell)
}
\arguments{
\item{bio02V}{Vector or single-column matrix of bio02 values.}

\item{bio07V}{Vector or single-column matrix of bio07 values.}

\item{cell}{Vector of original cell IDs.}
}
\value{
Matrix with "bio03", "cell".
}
\description{
Calculates (bio02 / bio07) * 100.
}
\keyword{internal}
