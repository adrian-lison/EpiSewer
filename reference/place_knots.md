# Places knots for fitting B-splines to a time series

Places knots for fitting B-splines to a time series

## Usage

``` r
place_knots(
  length_R,
  forecast_horizon,
  knot_distance,
  partial_window,
  partial_generation,
  fix_forecast = FALSE
)
```

## Arguments

- length_R:

  Length of the Rt time series

- forecast_horizon:

  Number of days to forecast

- knot_distance:

  Normal distance between knots

- partial_window:

  Window with only partial data towards the present in which the knot
  distances will be shorter. This is to avoid erroneous extrapolation of
  splines towards the present.

- partial_generation:

  A certain quantile of the generation time distribution. The first
  "fully informed" knot is placed at that distance after the partial
  window. The reason for this is that while the number of infections is
  already sufficiently informed after the partial window, the trend in
  infections is typically very volatile, so we place the first variable
  knot after approximately one generation.

## Value

A vector with knot positions.

## Details

The splines are forced to a zero slope towards the present and for
forecasting.
