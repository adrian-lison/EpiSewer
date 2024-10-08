% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_shedding.R
\name{load_per_case_assume}
\alias{load_per_case_assume}
\title{Assume the average load per case}
\usage{
load_per_case_assume(load_per_case = NULL, modeldata = modeldata_init())
}
\arguments{
\item{load_per_case}{The assumed average total shedding load per case. Must
have the same unit as the numerator of the concentration unit. For example,
if concentration is measured in gc/mL (gc = gene copies), then
\code{load_per_case} should also be in gc.}

\item{modeldata}{A \code{modeldata} object to which the above model specifications
should be added. Default is an empty model given by \code{\link[=modeldata_init]{modeldata_init()}}. Can
also be an already partly specified model returned by other \code{EpiSewer}
modeling functions.}
}
\value{
A \code{modeldata} object containing data and specifications of the model
to be fitted. Can be passed on to other \code{EpiSewer} modeling functions to
add further data and model specifications.

The \code{modeldata} object also includes information about parameter
initialization (\code{init}), meta data (\code{.metainfo}), and checks to be
performed before model fitting (\code{.checks}).
}
\description{
This option assumes an average total shedding load per case. In
the \code{EpiSewer} model, this serves as a scaling factor describing how
many pathogen particles are shed by the average infected individual overall
and how much of this is detectable at the sampling site. This depends
both on biological factors as well as on the specific sewage system.
}
\seealso{
{Helper for finding a suitable load per case assumption:}
\code{\link[=suggest_load_per_case]{suggest_load_per_case()}}

Other {load per case functions}: 
\code{\link{load_per_case_calibrate}()}
}
\concept{{load per case functions}}
