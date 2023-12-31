% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_sampling.R
\name{sample_effects_estimate_matrix}
\alias{sample_effects_estimate_matrix}
\title{Estimate sample effects using a design matrix}
\usage{
sample_effects_estimate_matrix(
  design_matrix,
  effect_prior_mu = 0,
  effect_prior_sigma = 1,
  modeldata = modeldata_init()
)
}
\arguments{
\item{design_matrix}{A design matrix with different covariates potentially
influencing sample concentration. The matrix must have one row for each
modeled day. See \code{\link[stats:model.matrix]{stats::model.matrix()}} for construction of design
matrices.}

\item{effect_prior_mu}{Prior (mean) on the regression coefficients.}

\item{effect_prior_sigma}{Prior (standard deviation) on the regression
coefficients.}

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
This option uses a log-linear regression model to estimate
effects of sample covariates on the concentration. Concentrations can be
influenced by sampling-related external factors, for example the time
between sampling and shipping to the lab (age-of-sample effect), or
different sampling or storage methods.
}
\details{
\code{EpiSewer} will fit a fixed-effects log-linear model, random effects
are currently not supported.

The priors of this component have the following functional form:
\itemize{
\item regression coefficients: \code{Normal}
}
}
\seealso{
Other {sample effect models}: 
\code{\link{sample_effects_estimate_weekday}()},
\code{\link{sample_effects_none}()}
}
\concept{{sample effect models}}
