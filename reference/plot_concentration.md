# Plot predicted concentration

Plots measured and/or predicted concentrations over time. The predicted
concentrations are taken from a fitted `EpiSewer` model.

## Usage

``` r
plot_concentration(
  results = NULL,
  measurements = NULL,
  flows = NULL,
  include_noise = TRUE,
  normalized = FALSE,
  median = FALSE,
  mark_outliers = FALSE,
  forecast = TRUE,
  forecast_horizon = NULL,
  date_margin_left = 0,
  date_margin_right = 1,
  facet_models = FALSE,
  facet_direction = "rows",
  base_model = "",
  model_levels = NULL,
  concentration_col = "concentration",
  flow_col = "flow",
  date_col = "date",
  outlier_col = "is_outlier",
  type = "time",
  intervals = c(0.5, 0.95),
  obs_size = 1.5,
  obs_shape = 4,
  obs_forecast_shape = 8
)
```

## Arguments

- results:

  Results object returned by
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  after model fitting. Can also be a named `list` with results from
  different model runs, in which case all results are plotted together
  and distinguished by colors.

- measurements:

  A `data.frame` with observed measurements, which will be plotted
  alongside the predicted values (useful to assess model fit). Can also
  be used with `results = NULL`, in which case only the observed
  measurements are plotted.

- flows:

  A `data.frame` with observed flows. If `normalized = TRUE`, this will
  be used to normalize the estimated and observed concentrations to the
  median flow.

- include_noise:

  If `TRUE` (default), concentrations including measurement noise are
  shown. If `FALSE`, only the expected concentrations are shown.

- normalized:

  If `TRUE`, the estimated and observed concentrations are normalized to
  the median flow. Thereby, the noise in concentrations that is due to
  variation in flow is essentially removed. This is especially useful
  for assessing forecasting performance when the future flow values were
  not known during model fitting.

- median:

  Should the estimated median be shown, or only the credible intervals?
  Default is `FALSE` to avoid over-interpretation of the median.

- mark_outliers:

  If `TRUE`, outliers in the `measurements` are highlighted in red. See
  also argument `outlier_col` below. Currently, this only works with
  manually labeled outliers. Outliers detected by the model can be found
  in `summary$outliers`.

- forecast:

  Should forecasted concentrations be shown? Default is true. This
  requires that the model was fitted with a forecast horizon, see
  [`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md).
  Because concentrations depend on the flow volume, accurate
  concentration forecasts are only possible if flow values beyond the
  estimation date are provided. This is typically only possible when
  estimating retrospectively. If flow values beyond the estimation date
  are not available, the forecasted concentrations will be based on the
  median flow volume.

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
  be negative to cut off some of the latest dates. By default, this is 1
  to improve the visibility of the latest date.

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

- concentration_col:

  Name of the column in the measurements `data.frame` which contains the
  measured concentration that should be plotted.

- date_col:

  Name of the date column in the measurements `data.frame`.

- outlier_col:

  Name of a logical column in the measurements `data.frame` which
  identifies outlier measurements (for example added by
  [`mark_outlier_spikes_median()`](https://adrian-lison.github.io/EpiSewer/reference/mark_outlier_spikes_median.md)).

- type:

  If "time" (default), the concentration time series is plotted. If
  "pp_check", the observations are ordered by concentration and plotted
  against the predicted concentration (useful for posterior predictive
  checks).

- obs_size:

  Size of the observed concentration points. Default is 1.5.

- obs_shape:

  Shape of the observed concentration points. Default is 4.

- obs_forecast_shape:

  Shape of the observed concentration points for forecasted values.
  Default is 8.

## Value

A `ggplot` object showing predicted and observed concentrations over
time. Can be further manipulated using `ggplot2` functions to adjust
themes and scales, and to add further geoms.

## Details

When plotting a posterior predictive check (`type="pp_check"`), each
observed concentration is shown together with the 95% Credible Interval
of the posterior predictive distribution for that observation. The
observations are ordered by concentration, which helps to visualize bias
and variance of the predictions as a function of the concentration.
