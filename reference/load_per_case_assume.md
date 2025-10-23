# Assume the average load per case

This option assumes an average total shedding load per case. In the
`EpiSewer` model, this serves as a scaling factor describing how many
pathogen particles are shed by the average infected individual overall
and how much of this is detectable at the sampling site. This depends
both on biological factors as well as on the specific sewage system.

## Usage

``` r
load_per_case_assume(load_per_case = NULL, modeldata = modeldata_init())
```

## Arguments

- load_per_case:

  The assumed average total shedding load per case. Must have the same
  unit as the numerator of the concentration unit. For example, if
  concentration is measured in gc/mL (gc = gene copies), then
  `load_per_case` should also be in gc.

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

Helper for finding a suitable load per case assumption:
[`suggest_load_per_case()`](https://adrian-lison.github.io/EpiSewer/reference/suggest_load_per_case.md)

Other load per case functions:
[`load_per_case_calibrate()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_calibrate.md)
