# Specify observation data

Specify wastewater observation data such as concentration measurements
and flows. This is a convenience function to collect all observation
data in one object.

## Usage

``` r
sewer_data(measurements = NULL, flows = NULL, cases = NULL, ...)
```

## Arguments

- measurements:

  A `data.frame` with measured concentrations of the pathogen of
  interest. Will be automatically passed to
  [`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md).

- flows:

  A `data.frame` with wastewater flow volumes at the sampling site for
  each day. Will be automatically passed to
  [`flows_observe()`](https://adrian-lison.github.io/EpiSewer/reference/flows_observe.md).

- cases:

  A `data.frame` of case numbers with each row representing one day.
  Must have at least a column with dates and a column with case numbers.
  This data is not used for model fitting, but for calibration of the
  `load_per_case` assumption. Will be automatically passed to
  [`load_per_case_calibrate()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_calibrate.md).

- ...:

  Further observations to be passed to
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md).

## Value

A `list` with all observations supplied. Can be passed to the `data`
argument in
[`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md).
