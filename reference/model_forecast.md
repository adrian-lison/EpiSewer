# Forecasting module

This module function is used to specify the components of the `forecast`
module in `EpiSewer`.

Each component can be specified using one or several helper functions
(see available options below). See the documentation of the individual
helper functions to adjust various settings.

## Usage

``` r
model_forecast(
  horizon = horizon_none(),
  damping = damping_assume(damping = 0.95)
)
```

## Arguments

- horizon:

  The forecast horizon. How many days into the future should EpiSewer
  forecast? Note that this functionality is intended for short-term
  forecasts. Projections over longer horizons can be highly inaccurate.
  Available options:

  - [`horizon_none()`](https://adrian-lison.github.io/EpiSewer/reference/horizon_none.md)

  - [`horizon_assume()`](https://adrian-lison.github.io/EpiSewer/reference/horizon_assume.md)

- damping:

  EpiSewer dampens the forecast of Rt so that trends in transmission
  will level off after some time. This prevents unrealistic
  extrapolation of transmission dynamics. Available options:

  - [`damping_none()`](https://adrian-lison.github.io/EpiSewer/reference/damping_none.md)

  - [`damping_assume()`](https://adrian-lison.github.io/EpiSewer/reference/damping_assume.md)

## Value

A `modeldata` object containing the data and specifications of the
`forecast` module.

## Details

Forecasts account for the estimated variation of transmission dynamics
over time and therefore tend to become more uncertain at longer forecast
horizons. However, it is important to keep in mind that depending on the
Rt model used, EpiSewer will project the current transmission dynamics
to continue unchanged (when using
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md),
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md),
[`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md))
or according to a (dampened) linear trend (when using
[`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md)
and
[`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md)).
This assumption can be violated by various factors such as depletion of
susceptible individuals, changes in behavior, or public health
interventions.

## See also

Other module functions:
[`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md),
[`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md),
[`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md),
[`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md),
[`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md)
