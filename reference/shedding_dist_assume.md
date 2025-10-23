# Assume a shedding load distribution

This option assumes a fixed shedding load distribution for the shedding
model in `EpiSewer`.

## Usage

``` r
shedding_dist_assume(
  shedding_dist = NULL,
  shedding_reference = NULL,
  modeldata = modeldata_init()
)
```

## Arguments

- shedding_dist:

  A numeric vector representing a discrete shedding load distribution,
  with elements describing the share of load shed 0 days, 1 day, 2 days,
  and so on after the start of shedding.

- shedding_reference:

  Is the shedding load distribution relative to the day of `"infection"`
  or the day of `"symptom_onset"`? This is important because shedding
  load distributions provided in the literature are sometimes by days
  since infection and sometimes by days since symptom onset. If
  `shedding_reference="symptom_onset"`, EpiSewer also needs information
  about the incubation period distribution (see
  [`incubation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/incubation_dist_assume.md)).

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

Helpers to discretize continuous probability distributions:
[`get_discrete_gamma()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_gamma.md),
[`get_discrete_lognormal()`](https://adrian-lison.github.io/EpiSewer/reference/get_discrete_lognormal.md)
