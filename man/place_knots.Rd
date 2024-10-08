% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_helpers.R
\name{place_knots}
\alias{place_knots}
\title{Places knots for fitting B-splines to a time series}
\usage{
place_knots(ts_length, knot_distance, partial_window = 30)
}
\arguments{
\item{ts_length}{Length of the time series}

\item{knot_distance}{Normal distance between knots}

\item{partial_window}{Window with only partial data towards the present in
which the knot distances will be shorter. This is to avoid erroneous
extrapolation of splines towards the present.}
}
\value{
A vector with knot positions.
}
\description{
Places knots for fitting B-splines to a time series
}
\keyword{internal}
