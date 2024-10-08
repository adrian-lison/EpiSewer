% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_infections.R
\name{model_infections}
\alias{model_infections}
\title{Model the infection process}
\usage{
model_infections(
  generation_dist = generation_dist_assume(),
  R = R_estimate_splines(),
  seeding = seeding_estimate_rw(),
  infection_noise = infection_noise_estimate()
)
}
\arguments{
\item{generation_dist}{Generation time distribution. The intrinsic
distribution of the time between infection of a primary case and infection
of its secondary cases. Modeling options:
\itemize{
\item \code{\link[=generation_dist_assume]{generation_dist_assume()}}
}}

\item{R}{Effective reproduction number over time. This is the main parameter
of interest estimated by \code{EpiSewer}. \code{R} is smoothed using a time series
smoothing prior. Currently supported are: random walk (rw), exponential
smoothing (ets), and smoothing splines. Modeling options:
\itemize{
\item \code{\link[=R_estimate_splines]{R_estimate_splines()}}
\item \code{\link[=R_estimate_ets]{R_estimate_ets()}}
\item \code{\link[=R_estimate_approx]{R_estimate_approx()}}
\item \code{\link[=R_estimate_rw]{R_estimate_rw()}}
}}

\item{seeding}{Seeding of initial infections. The renewal model used by
\code{EpiSewer} requires a seeding phase of the length of the maximum generation
time. For these initial infections, a simple seeding model instead of the
renewal model must be used. Modeling options:
\itemize{
\item \code{\link[=seeding_estimate_constant]{seeding_estimate_constant()}}
\item \code{\link[=seeding_estimate_rw]{seeding_estimate_rw()}}
}}

\item{infection_noise}{Noise in the infection process. \code{EpiSewer} implements
a stochastic infection model, i.e. allows for variation in the number of new
infections generated at each time step. This accounts for stochastic
uncertainty in the infection process and often speeds up model fitting.
Modeling options:
\itemize{
\item \code{\link[=infection_noise_none]{infection_noise_none()}}
\item \code{\link[=infection_noise_estimate]{infection_noise_estimate()}}
}}
}
\value{
A \code{modeldata} object containing the data and specifications of the
\code{infections} module.
}
\description{
This module function is used to specify the components of the
\code{infections} module in \code{EpiSewer}.

Each component can be specified using one or several helper
functions (see available options below). See the documentation of the
individual helper functions to adjust model priors and further settings.
}
\seealso{
Other {module functions}: 
\code{\link{model_forecast}()},
\code{\link{model_measurements}()},
\code{\link{model_sampling}()},
\code{\link{model_sewage}()},
\code{\link{model_shedding}()}
}
\concept{{module functions}}
