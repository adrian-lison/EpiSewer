% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_measurements.R
\name{LOD_estimate_ddPCR}
\alias{LOD_estimate_ddPCR}
\title{Estimate a limit of detection model for digital droplet PCR data}
\usage{
LOD_estimate_ddPCR(modeldata = modeldata_init())
}
\arguments{
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
Pathogen concentrations below a certain threshold may not be
detectable and thus erroneously measured as 0. This option adjusts for a
limit of detection based on the statistical properties of digital droplet
PCR (ddPCR) and includes zero measurements in the likelihood.

In effect, zero measurements provide a signal that the
concentration in the respective sample was likely below the limit of
detection, but we don't know what the exact concentration was.
}
\details{
The limit of detection is modeled using a hurdle model. The model
uses the number of droplets in the ddPCR reaction and the concentration
scaling factor as defined and estimated by \code{\link[=noise_estimate_ddPCR]{noise_estimate_ddPCR()}}. It can
therefore only be used together with \code{noise = noise_estimate_ddPCR()} in
\code{\link[=model_measurements]{model_measurements()}}.
}
\seealso{
Other {LOD models}: 
\code{\link{LOD_assume}()},
\code{\link{LOD_none}()}
}
\concept{{LOD models}}
