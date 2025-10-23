# Compute crude estimate of load by date of infections

The function smoothes the measured loads using a loess smoother on the
log scale and shifts them backward by the mean total delay to obtain a
crude estimate of the loads by date of infection. This is useful for
different calibration purposes, such as determining a suitable
`load_per_case` assumption.

## Usage

``` r
get_load_curve_crude(
  measured_concentrations,
  measure_to_sample,
  sample_to_date,
  flow,
  total_delay_dist,
  max_shift = NULL,
  T_start_date,
  impute_zero = NA,
  impute_zero_runs = FALSE,
  interpolate = TRUE,
  loess_window = 56,
  plot_smoothed_curve = FALSE
)
```

## Arguments

- max_shift:

  How far into the past should the time series be expanded? If NULL
  (default), only expands by the mean total delay. If `max_shift` is
  supplied and larger than he mean total delay, LOESS extrapolation is
  used.

- impute_zero:

  How should zero measurements be imputed? If NA (default), they are
  removed. This is the imputation value for the measured concentration.
  Imputed loads are calculated as `impute_zero * flow`.

- impute_zero_runs:

  Should runs of zero measurements be imputed by accounting for how many
  measurements in a row have been zero? This is based on the posterior
  expectation of the concentration for a number of zero measurements. In
  effect, the imputed concentration will get smaller and smaller during
  a run of zeros.

- interpolate:

  Should a log-linear interpolation of missing dates be performed? This
  makes the loess smoothing window more consistent.

- loess_window:

  Length of the loess smoothing window in days. Will be converted to the
  span parameter.

- plot_smoothed_curve:

  Should the smoothed curve be plotted for inspection purposes?

## Value

A `data.frame` with columns `date_index`, `date`, and `load`.
