# Plot the effective reproduction number

Plots the effective reproduction number over time from one or several
fitted `EpiSewer` models.

## Usage

``` r
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
  intervals = c(0.5, 0.95),
  color_by_chain = FALSE
)
```

## Arguments

- results:

  Results object returned by
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  after model fitting. Can also be a named `list` with results from
  different model runs, in which case all results are plotted together
  and distinguished by colors.

- draws:

  If `FALSE` (default), 50% and 95% Bayesian credible intervals are
  shown. If `TRUE`, exemplary posterior samples are shown in a
  "spaghetti plot" style.

- ndraws:

  Number of different samples to show if `draws=TRUE`.

- median:

  Should the estimated median be shown, or only the credible intervals?
  Default is `FALSE` to avoid over-interpretation of the median.

- seeding:

  Should Rt from the seeding phase be shown as well? Default is `FALSE`.

- forecast:

  Should forecasted Rt values be shown? Default is true. This requires
  that the model was fitted with a forecast horizon, see
  [`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md).

- forecast_horizon:

  How many days into the future should forecasts be plotted? Note that
  this is restricted by the forecast horizon specified during model
  fitting, see
  [`horizon_assume()`](https://adrian-lison.github.io/EpiSewer/reference/horizon_assume.md).

- date_margin_left:

  By how many days into the past should the plot be expanded? Can also
  be negative to cut off some of the earliest dates.

- date_margin_right:

  By how many days into the future should the plot be expanded? Can also
  be negative to cut off some of the latest dates.

- facet_models:

  Should the plot be faceted by model? Default is `FALSE`.

- facet_direction:

  How should the facetting be done? Either in different "rows" (default)
  or in different "columns".

- base_model:

  Name of the base model (in the named list provided to `results`) which
  should be compared to the other models. This model will be plotted in
  black and will not be part of the legend.

- model_levels:

  A `character` vector with the names of the models to be included. The
  colors and legend will be ordered according to the order in
  `model_levels`.

- color_by_chain:

  If `TRUE`, and `draws=TRUE`, individual samples are colored by chain.
  This is useful for checking mixing of the chains with respect to the
  estimated trajectories. Only available when plotting a single model.

## Value

A `ggplot` object showing the time series of estimated Rt, either with
credible intervals or as "spaghetti plot". Can be further manipulated
using `ggplot2` functions to adjust themes and scales, and to add
further geoms.
