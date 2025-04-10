% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_infections.R
\name{seeding_estimate_rw}
\alias{seeding_estimate_rw}
\title{Estimate seeding infections using a random walk model}
\usage{
seeding_estimate_rw(
  intercept_prior_q5 = NULL,
  intercept_prior_q95 = NULL,
  rel_change_prior_mu = 0.05,
  rel_change_prior_sigma = 0.025,
  extend = TRUE,
  modeldata = modeldata_init()
)
}
\arguments{
\item{intercept_prior_q5}{Prior (5\% quantile) on the initial number of
infections. Can be interpreted as an approximate lower bound. If NULL
(default), this is computed from a crude empirical estimate of the number
of cases (see details).}

\item{intercept_prior_q95}{Prior (95\% quantile) on the initial number of
infections. Can be interpreted as an approximate upper bound. If NULL
(default), this is computed from a crude empirical estimate of the number
of cases (see details).}

\item{rel_change_prior_mu}{Prior (mean) on the relative change rate of the
geometric random walk during the seeding phase. The default value (0.05)
assumes that daily changes are +-5\% on expectation and likely less than
+-10\% per day.}

\item{rel_change_prior_sigma}{Prior (standard deviation) on the relative
change rate of the geometric random walk during the seeding phase. This
expresses your uncertainty about the change rate. The default value (0.025)
assumes that the daily change rate could be 5\% points higher or lower than
your prior mean. For example, if \code{rel_change_prior_mu = 0.05} and
\code{rel_change_prior_sigma = 0.025}, this means you expect the daily change
rate to be between 0 (0\%) and 0.1 (10\%).}

\item{extend}{Should the seeding phase be extended when concentrations are
very low at the start of the measurement time series? If \code{TRUE}, then the
seeding phase will be extended to the first date with three consecutive
detects (i.e. non-zero measurements). The reproduction number will only be
modeled from that date onward. This option often makes sense, as infection
numbers are typically very low during a period with many non-detects, which
can lead to sampling problems when estimating Rt. If you nevertheless want
Rt estimates also for this period, you can use \code{extend = FALSE}. Note
though that estimated reproduction numbers are not necessarily meaningful
during periods with very low infection numbers, as transmission dynamics
may be dominated by chance events and importations.}

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
a geometric random walk to model these seeding infections.
}
\details{
The seeding phase has the length of the maximum generation time
(during this time, the renewal model cannot be applied). Traditionally,
seeding refers to the first few (potentially imported) infections of an
epidemic, but depending on what time period the model is fitted to, this
may also cover a different phase with stronger growth dynamics.

If \code{intercept_prior_q5} or \code{intercept_prior_q95} are not specified
by the user, \code{EpiSewer} will compute a rough median empirical estimate of
the number of cases using the supplied wastewater measurements and shedding
assumptions, and then infer the missing quantiles based on this. If none of
the quantiles are provided, they are set to be roughly 1/10 and 10 times
the empirical median estimate. We note that this is a violation of Bayesian
principles (data must not be used to inform priors) - but a neglectable
one, since it only ensures that the seeding is modeled on the right order
of magnitude and does not have relevant impacts on later Rt estimates.

The priors of this component have the following functional form:
\itemize{
\item intercept of the random walk (log scale): \code{Normal}
\item standard deviation of the random walk (log scale): \verb{Truncated normal}
The priors for these parameters are determined based on the user-supplied
arguments, using appropriate transformations and the two-sigma-rule of
thumb.
}
}
