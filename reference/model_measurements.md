# Model the measurement process

This module function is used to specify the components of the
`measurements` module in `EpiSewer`.

Each component can be specified using one or several helper functions
(see available options below). See the documentation of the individual
helper functions to adjust model priors and further settings.

## Usage

``` r
model_measurements(
  concentrations = concentrations_observe(),
  noise = noise_estimate(),
  LOD = LOD_none()
)
```

## Arguments

- concentrations:

  Concentration measurements from wastewater samples. Modeling options:

  - [`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md)

- noise:

  Measurement noise due to unexplained variation in sampling and lab
  analysis. Modeling options:

  - [`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md)

  - [`noise_estimate_constant_var()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_constant_var.md)

  - [`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)

- LOD:

  Limit of detection. Concentrations below a certain threshold may not
  be detectable and thus erroneously measured as 0. `EpiSewer` can
  adjust for the limit of detection using a hurdle model. Modeling
  options:

  - [`LOD_none()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_none.md)

  - [`LOD_assume()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_assume.md)

  - [`LOD_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)

## Value

A `modeldata` object containing the data and specifications of the
`measurements` module.

## See also

Other module functions:
[`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md),
[`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md),
[`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md),
[`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md),
[`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md)
