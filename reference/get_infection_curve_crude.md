# Compute crude estimate of infections by date of infections

The function smoothed the measured loads using a loess smoother on the
log scale and shifts them backward by the mean total delay to obtain a
crude estimate of the loads by date of infection. This is useful for
different calibration purposes, such as determining a suitable
`load_per_case` assumption.

## Usage

``` r
get_infection_curve_crude(load_curve_crude, load_per_case)
```

## Arguments

- load_curve_crude:

  A `data.frame` with columns `date_index`, `date`, and `load` as
  produced by
  [`get_load_curve_crude()`](https://adrian-lison.github.io/EpiSewer/reference/get_load_curve_crude.md).

- load_per_case:

  Load per case in gene copies.

## Value

A `data.frame` with columns `date_index`, `date`, and `load`.
