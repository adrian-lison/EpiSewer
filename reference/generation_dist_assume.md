# Assume a generation time distribution

This option assumes a fixed generation time distribution for the renewal
model in `EpiSewer`.

## Usage

``` r
generation_dist_assume(generation_dist = NULL, modeldata = modeldata_init())
```

## Arguments

- generation_dist:

  A numeric vector representing a discrete generation time distribution,
  starting with the probability for a generation time of 1 day, 2 days,
  3 days, and so on (a generation time of 0 days is excluded).

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

The generation time distribution here refers to the intrinsic
distribution of the time between infection of a primary case and
infection of its secondary cases. It is disease-specific and typically
obtained from literature.

## See also

Helpers to discretize continuous probability distributions:
[`get_discrete_gamma()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_gamma.md),
[`get_discrete_gamma_shifted()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_gamma_shifted.md),
[`get_discrete_lognormal()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_lognormal.md)
