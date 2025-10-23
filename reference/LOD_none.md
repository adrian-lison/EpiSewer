# Do not model a limit of detection

This option drops all zero measurements from the likelihood and does not
explicitly model a limit of detection.

## Usage

``` r
LOD_none(modeldata = modeldata_init())
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

## Details

Dropping zero measurements implicitly assumes that any level of
concentration can theoretically lead to a zero measurement. Since the
corresponding observations are dropped, this will discard information,
but not bias estimates.

## See also

Other LOD models:
[`LOD_assume()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_assume.md),
[`LOD_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
