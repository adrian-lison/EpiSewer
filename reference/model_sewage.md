# Model the sewage process

This module function is used to specify the components of the `sewage`
module in `EpiSewer`.

Each component can be specified using one or several helper functions
(see available options below). See the documentation of the individual
helper functions to adjust model priors and further settings.

## Usage

``` r
model_sewage(flows = flows_observe(), residence_dist = residence_dist_assume())
```

## Arguments

- flows:

  Daily flow volumes at the sampling site. The flow can change due to
  rainfall or industrial discharge, and directly influences pathogen
  concentrations in the wastewater. Modeling options:

  - [`flows_observe()`](https://adrian-lison.github.io/EpiSewer/reference/flows_observe.md)

  - [`flows_assume()`](https://adrian-lison.github.io/EpiSewer/reference/flows_assume.md)

- residence_dist:

  Sewer residence time distribution for pathogen particles. By default,
  `EpiSewer` assumes that particles arrive at the sampling site within
  the day of shedding. However, for larger sewage systems, particles may
  travel longer than a day depending on where and when they were shed
  into the wastewater. Modeling options:

  - [`residence_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/residence_dist_assume.md)

## Value

A `modeldata` object containing the data and specifications of the
`sewage` module.

## See also

Other module functions:
[`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md),
[`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md),
[`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md),
[`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md),
[`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md)
