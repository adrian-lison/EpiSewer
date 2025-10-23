# Plot the estimated load

Plots the estimated load in wastewater over time from a fitted
`EpiSewer` model.

## Usage

``` r
plot_load(
  results,
  median = FALSE,
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
```

## Arguments

- results:

  Results object returned by
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  after model fitting. Can also be a named `list` with results from
  different model runs, in which case all results are plotted together
  and distinguished by colors.

- median:

  Should the estimated median be shown, or only the credible intervals?
  Default is `FALSE` to avoid over-interpretation of the median.

- forecast:

  Should forecasted loads be shown? Default is true. This requires that
  the model was fitted with a forecast horizon, see
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

## Value

A `ggplot` object showing the estimated load over time. Can be further
manipulated using `ggplot2` functions to adjust themes and scales, and
to add further geoms.
