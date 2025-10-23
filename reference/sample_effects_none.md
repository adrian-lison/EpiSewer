# Do not model sample effects

This option does not model effects of sample covariates on the
concentrations.

## Usage

``` r
sample_effects_none(modeldata = modeldata_init())
```

## Arguments

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

## See also

Other sample effect models:
[`sample_effects_estimate_matrix()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_matrix.md),
[`sample_effects_estimate_weekday()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_weekday.md)
