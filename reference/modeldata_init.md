# Construct an unspecified EpiSewer model

Construct an unspecified EpiSewer model

## Usage

``` r
modeldata_init()
```

## Value

A `modeldata` object containing data and specifications of the model to
be fitted. Can be passed on to other `EpiSewer` modeling functions to
add further data and model specifications.

The `modeldata` object also includes information about parameter
initialization (`.init`), meta data (`.metainfo`), and checks to be
performed before model fitting (`.checks`).
