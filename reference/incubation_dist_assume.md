# Assume an incubation period distribution

This option assumes a fixed incubation period distribution for the
shedding model in `EpiSewer`.

## Usage

``` r
incubation_dist_assume(incubation_dist = NULL, modeldata = modeldata_init())
```

## Arguments

- incubation_dist:

  A numeric vector representing a discrete incubation period
  distribution, starting with the probability for an incubation period
  of 0 days, 1 day, 2 days, and so on.

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

The incubation period is the time between infection and symptom onset.
This assumption is used when `shedding_reference="symptom_onset"`, i.e.
to support shedding load distributions referenced by days since symptom
onset.

## See also

Helpers to discretize continuous probability distributions:
[`get_discrete_gamma()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_gamma.md),
[`get_discrete_lognormal()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_lognormal.md)
