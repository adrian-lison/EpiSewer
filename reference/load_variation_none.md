# Do not model individual-level load variation

This option assumes that there is no variation in the total load shed
per case, i.e. the individual shedding load is fixed to the average
shedding load.

## Usage

``` r
load_variation_none(modeldata = modeldata_init())
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

Other load variation models:
[`load_variation_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/load_variation_estimate.md)
