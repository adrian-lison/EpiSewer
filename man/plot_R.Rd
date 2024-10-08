% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sewer_plot.R
\name{plot_R}
\alias{plot_R}
\title{Plot the effective reproduction number}
\usage{
plot_R(
  results,
  draws = FALSE,
  ndraws = NULL,
  median = FALSE,
  seeding = FALSE,
  forecast = TRUE,
  forecast_horizon = NULL,
  date_margin_left = 0,
  date_margin_right = 0,
  facet_models = FALSE,
  facet_direction = "rows",
  base_model = "",
  model_levels = NULL,
  intervals = c(0.5, 0.95)
)
}
\arguments{
\item{results}{Results object returned by \code{\link[=EpiSewer]{EpiSewer()}} after model fitting.
Can also be a named \code{list} with results from different model runs, in which
case all results are plotted together and distinguished by colors.}

\item{draws}{If \code{FALSE} (default), 50\% and 95\% Bayesian credible intervals
are shown. If \code{TRUE}, exemplary posterior samples are shown in a "spaghetti
plot" style.}

\item{ndraws}{Number of different samples to show if \code{draws=TRUE}.}

\item{median}{Should the estimated median be shown, or only the credible
intervals? Default is \code{FALSE} to avoid over-interpretation of the median.}

\item{seeding}{Should Rt from the seeding phase be shown as well? Default is
\code{FALSE}.}

\item{forecast}{Should forecasted Rt values be shown? Default is true. This
requires that the model was fitted with a forecast horizon, see
\code{\link[=model_forecast]{model_forecast()}}.}

\item{forecast_horizon}{How many days into the future should forecasts be
plotted? Note that this is restricted by the forecast horizon specified
during model fitting, see \code{\link[=horizon_assume]{horizon_assume()}}.}

\item{date_margin_left}{By how many days into the past should the plot be
expanded? Can also be negative to cut off some of the earliest dates.}

\item{date_margin_right}{By how many days into the future should the plot be
expanded? Can also be negative to cut off some of the latest dates.}

\item{facet_models}{Should the plot be faceted by model? Default is \code{FALSE}.}

\item{facet_direction}{How should the facetting be done? Either in different
"rows" (default) or in different "columns".}

\item{base_model}{Name of the base model (in the named list provided to
\code{results}) which should be compared to the other models. This model will be
plotted in black and will not be part of the legend.}

\item{model_levels}{A \code{character} vector with the names of the models to be
included. The colors and legend will be ordered according to the order in
\code{model_levels}.}
}
\value{
A ggplot object showing the time series of estimated Rt, either with
credible intervals or as "spaghetti plot". Can be further manipulated using
\link{ggplot2} functions to adjust themes and scales, and to add further geoms.
}
\description{
Plots the effective reproduction number over time from one or
several fitted \code{EpiSewer} models.
}
