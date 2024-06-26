% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sewer_plot.R
\name{plot_load}
\alias{plot_load}
\title{Plot the estimated load}
\usage{
plot_load(
  results,
  median = FALSE,
  date_margin_left = 0,
  date_margin_right = 0,
  facet_models = FALSE,
  base_model = "",
  model_levels = NULL
)
}
\arguments{
\item{results}{Results object returned by \code{\link[=EpiSewer]{EpiSewer()}} after model fitting.
Can also be a named \code{list} with results from different model runs, in which
case all results are plotted together and distinguished by colors.}

\item{median}{Should the estimated median be shown, or only the credible
intervals? Default is \code{FALSE} to avoid over-interpretation of the median.}

\item{date_margin_left}{By how many days into the past should the plot be
expanded? Can also be negative to cut off some of the earliest dates.}

\item{date_margin_right}{By how many days into the future should the plot be
expanded? Can also be negative to cut off some of the latest dates.}

\item{facet_models}{Should the plot be faceted by model? Default is \code{FALSE}.}

\item{base_model}{Name of the base model (in the named list provided to
\code{results}) which should be compared to the other models. This model will be
plotted in black and will not be part of the legend.}

\item{model_levels}{A \code{character} vector with the names of the models to be
included. The colors and legend will be ordered according to the order in
\code{model_levels}.}
}
\value{
A ggplot object showing the estimated load over time. Can be further
manipulated using \link{ggplot2} functions to adjust themes and scales, and to
add further geoms.
}
\description{
Plots the estimated load in wastewater over time from a fitted
\code{EpiSewer} model.
}
