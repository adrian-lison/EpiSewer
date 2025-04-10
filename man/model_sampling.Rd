% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_sampling.R
\name{model_sampling}
\alias{model_sampling}
\title{Model the sampling process}
\usage{
model_sampling(
  outliers = outliers_none(),
  sample_effects = sample_effects_none()
)
}
\arguments{
\item{outliers}{Outliers in concentrations. \code{EpiSewer} can automatically
identify independent spikes in the concentration time series and model them
as outliers to reduce the impact on transmission dynamic estimates.
Modeling options:
\itemize{
\item \code{\link[=outliers_none]{outliers_none()}}
\item \code{\link[=outliers_estimate]{outliers_estimate()}}
}}

\item{sample_effects}{Sample (batch) effects. The pathogen concentration in a
sample may be influenced by sampling-related external factors, for example
the time between sampling and shipping to the lab (age-of-sample effect),
or different sampling or storage methods. \code{EpiSewer} allows to estimate
such effects using covariates that describe differences between the
samples. Modeling options:
\itemize{
\item \code{\link[=sample_effects_none]{sample_effects_none()}}
\item \code{\link[=sample_effects_estimate_weekday]{sample_effects_estimate_weekday()}}
\item \code{\link[=sample_effects_estimate_matrix]{sample_effects_estimate_matrix()}}
}}
}
\value{
A \code{modeldata} object containing the data and specifications of the
\code{sampling} module.
}
\description{
This module function is used to specify the components of the
\code{sampling} module in \code{EpiSewer}.

Each component can be specified using one or several helper
functions (see available options below). See the documentation of the
individual helper functions to adjust model priors and further settings.
}
\seealso{
Other {module functions}: 
\code{\link{model_forecast}()},
\code{\link{model_infections}()},
\code{\link{model_measurements}()},
\code{\link{model_sewage}()},
\code{\link{model_shedding}()}
}
\concept{{module functions}}
