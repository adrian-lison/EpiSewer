% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_measurements.R
\name{noise_estimate_}
\alias{noise_estimate_}
\title{Estimate measurement noise (internal helper function)}
\usage{
noise_estimate_(
  replicates = FALSE,
  cv_prior_mu = 0,
  cv_prior_sigma = 1,
  cv_type = "constant",
  ddPCR_prior_droplets_mu = NULL,
  ddPCR_prior_droplets_sigma = NULL,
  ddPCR_droplets_fixed = NULL,
  ddPCR_prior_scaling_mu = NULL,
  ddPCR_prior_scaling_sigma = NULL,
  ddPCR_scaling_fixed = NULL,
  pre_replicate_cv_prior_mu = 0,
  pre_replicate_cv_prior_sigma = 1,
  modeldata = modeldata_init()
)
}
\arguments{
\item{replicates}{Should replicates be used to explicitly model variation
before the replication stage?}

\item{cv_prior_mu}{Prior (mean) on the coefficient of variation of
concentration measurements.}

\item{cv_prior_sigma}{Prior (standard deviation) on the coefficient of
variation of concentration measurements.}

\item{cv_type}{One out of "constant" (default), or "ddPCR". If ddPCR, the
coefficient of variation is modeled as a function of the expected
concentration according to the statistical properties of ddPCR. In
particular, this model predicts a higher coefficient of variation at
smaller concentrations, which often leads to a better model fit.}

\item{ddPCR_prior_droplets_mu}{Prior (mean) on the number of droplets in the
ddPCR reaction.}

\item{ddPCR_prior_droplets_sigma}{Prior (standard deviation) on the number of
droplets in the ddPCR reaction.}

\item{ddPCR_droplets_fixed}{If TRUE (default), the number of droplets is
fixed to the prior mean and not estimated. This is recommended if no
replicates are available.}

\item{ddPCR_prior_scaling_mu}{Prior (mean) on the concentration scaling
factor for the ddPCR reaction. The concentration scaling factor is the
droplet volume, scaled by the dilution of the wastewater in the ddPCR
reaction. See details for further explanation.}

\item{ddPCR_prior_scaling_sigma}{Prior (standard deviation) on the
concentration scaling factor for the ddPCR reaction.}

\item{ddPCR_scaling_fixed}{If TRUE, the concentration scaling factor is fixed
to the prior mean and not estimated.}

\item{pre_replicate_cv_prior_mu}{Prior (mean) on the coefficient of variation
of concentrations before the replication stage.}

\item{pre_replicate_cv_prior_sigma}{Prior (standard deviation) on the
coefficient of variation of concentrations before the replication stage.}

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
This option estimates the unexplained variation in wastewater
measurements. If multiple measurements (replicates) per sample are
provided, \code{EpiSewer} can also explicitly model variation before the
replication stage.

This helper function is called from \code{\link[=noise_estimate]{noise_estimate()}} and
\code{\link[=noise_estimate_ddPCR]{noise_estimate_ddPCR()}}. \code{\link[=noise_estimate]{noise_estimate()}} is a constant coefficient of
variation model, \code{\link[=noise_estimate_ddPCR]{noise_estimate_ddPCR()}} is a noise model specialized for
digital droplet PCR (\code{cv_type = "ddPCR"}), which may however
also work with other quantification methods such as qPCR.
}
