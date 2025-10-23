# Model the sampling process

This module function is used to specify the components of the `sampling`
module in `EpiSewer`.

Each component can be specified using one or several helper functions
(see available options below). See the documentation of the individual
helper functions to adjust model priors and further settings.

## Usage

``` r
model_sampling(
  outliers = outliers_none(),
  sample_effects = sample_effects_none()
)
```

## Arguments

- outliers:

  Outliers in concentrations. `EpiSewer` can automatically identify
  independent spikes in the concentration time series and model them as
  outliers to reduce the impact on transmission dynamic estimates.
  Modeling options:

  - [`outliers_none()`](https://adrian-lison.github.io/EpiSewer/reference/outliers_none.md)

  - [`outliers_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/outliers_estimate.md)

- sample_effects:

  Sample (batch) effects. The pathogen concentration in a sample may be
  influenced by sampling-related external factors, for example the time
  between sampling and shipping to the lab (age-of-sample effect), or
  different sampling or storage methods. `EpiSewer` allows to estimate
  such effects using covariates that describe differences between the
  samples. Modeling options:

  - [`sample_effects_none()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_none.md)

  - [`sample_effects_estimate_matrix()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_matrix.md)

  - [`sample_effects_estimate_weekday()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_weekday.md)

## Value

A `modeldata` object containing the data and specifications of the
`sampling` module.

## See also

Other module functions:
[`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md),
[`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md),
[`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md),
[`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md),
[`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md)
