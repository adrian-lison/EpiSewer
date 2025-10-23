# Specify the forecast horizon

This option specifies a fixed forecast horizon in days. EpiSewer will
produce concentration predictions and forecasts of all latent variables
such as Rt, infections and load until the end of the forecast horizon.

## Usage

``` r
horizon_assume(horizon, modeldata = modeldata_init())
```

## Arguments

- horizon:

  The forecast horizon in days. If 0, no forecasts are produced. Note
  that this functionality is intended for short-term forecasts.
  Projections over longer horizons can be highly inaccurate.

- modeldata:

  A `modeldata` object to which the above model specifications should be
  added. Default is an empty model given by
  [`modeldata_init()`](https://adrian-lison.github.io/EpiSewer/reference/modeldata_init.md).
  Can also be an already partly specified model returned by other
  `EpiSewer` modeling functions.

## Value

A `modeldata` object containing data and specifications of the model to
be fitted. Can be passed on to other `EpiSewer` modeling functions to
add further data and model specifications.

The `modeldata` object also includes information about parameter
initialization (`.init`), meta data (`.metainfo`), and checks to be
performed before model fitting (`.checks`).

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
