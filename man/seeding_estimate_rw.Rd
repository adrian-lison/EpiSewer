% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_infections.R
\name{seeding_estimate_rw}
\alias{seeding_estimate_rw}
\title{Estimate seeding infections}
\usage{
seeding_estimate_rw(
  intercept_prior_mu = NULL,
  intercept_prior_sigma = 1,
  sd_prior_mu = 0.05,
  sd_prior_sigma = 0.025,
  modeldata = modeldata_init()
)
}
\arguments{
\item{intercept_prior_mu}{Prior (mean) on the intercept of the seeding random
walk. If NULL (default), this is set to a crude empirical estimate of the
number of cases (see details).}

\item{intercept_prior_sigma}{Prior (standard deviation) on the intercept of
the seeding random walk.}

\item{sd_prior_mu}{Prior (mean) on the standard deviation of the seeding
random walk.}

\item{sd_prior_sigma}{Prior (standard deviation) on the standard deviation of
the seeding random walk.}

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
This option estimates initial infections at the start of the
modeled time period, when the renewal model cannot be applied yet. It uses
a random walk to model these seeding infections.
}
\details{
The seeding phase has the length of the maximum generation time
(during this time, the renewal model cannot be applied). Traditionally,
seeding refers to the first few (potentially imported) infections of an
epidemic, but depending on what time period the model is fitted to, this
may also cover a different phase with higher incidence levels.

If \code{intercept_prior_mu} is not specified by the user, \code{EpiSewer}
will set it to a rough estimate of the number of cases using the supplied
wastewater measurements and shedding assumptions. We note that this is a
violation of Bayesian principles (data must not be used to inform priors) -
but a neglectable one, since it only ensures that the seeding is modeled on
the right order of magnitude and does not have relevant impacts on later Rt
estimates.

The priors of this component have the following functional form:
\itemize{
\item intercept of the seeding random walk: \code{Normal}
\item standard deviation of the seeding random walk: \verb{Truncated normal}
}
}
