# Model the shedding process

This module function is used to specify the components of the `shedding`
module in `EpiSewer`.

Each component can be specified using one or several helper functions
(see available options below). See the documentation of the individual
helper functions to adjust model priors and further settings.

## Usage

``` r
model_shedding(
  incubation_dist = incubation_dist_assume(),
  shedding_dist = shedding_dist_assume(),
  load_per_case = load_per_case_calibrate(),
  load_variation = load_variation_estimate()
)
```

## Arguments

- incubation_dist:

  Incubation period distribution. `EpiSewer` uses this as a proxy for
  the time between infection and the start of shedding, as shedding load
  distributions in the literature are often given from symptom onset
  onwards. If the assumed shedding load distribution instead starts from
  the time of infection, the incubation period should be fixed to 0
  days. Modeling options:

  - [`incubation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/incubation_dist_assume.md)

- shedding_dist:

  Shedding load distribution. Describes how the total load shed by an
  individual is distributed over time (and therefore sums to 1).
  Modeling options:

  - [`shedding_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/shedding_dist_assume.md)

  - [`shedding_dist_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/shedding_dist_estimate.md)

- load_per_case:

  Average total load per case. This is a scaling factor that describes
  how many pathogen particles are shed by the average infected
  individual overall and how much of this is detectable at the sampling
  site. It depends both on biological factors as well as on the specific
  sewage system. Modeling options:

  - [`load_per_case_assume()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_assume.md)

- load_variation:

  Individual-level shedding load variation. The strength of shedding may
  vary between individuals. Modeling this variation can better account
  for uncertainty especially at low incidence. Modeling options:

  - [`load_variation_none()`](https://adrian-lison.github.io/EpiSewer/reference/load_variation_none.md)

  - [`load_variation_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/load_variation_estimate.md)

## Value

A `modeldata` object containing the data and specifications of the
`shedding` module.

## See also

Other module functions:
[`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md),
[`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md),
[`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md),
[`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md),
[`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md)
