# Estimate epidemiological parameters from wastewater measurements

`EpiSewer` estimates the effective reproduction number Rt and other
parameters of interest from wastewater measurements over time. This
function combines data, assumptions and modeling details to fit a
Bayesian hierarchical model to the measurements. The resulting posterior
estimates of Rt and other parameters can be plotted and further
analyzed.

The inputs to this function can be specified using various helper
functions (see below). It is best to define the inputs in advance as
separate variables and then pass them to `EpiSewer`.

## Usage

``` r
EpiSewer(
  data = sewer_data(),
  assumptions = sewer_assumptions(),
  measurements = model_measurements(),
  sampling = model_sampling(),
  sewage = model_sewage(),
  shedding = model_shedding(),
  infections = model_infections(),
  forecast = model_forecast(),
  fit_opts = set_fit_opts(),
  results_opts = set_results_opts(),
  run_fit = TRUE
)
```

## Arguments

- data:

  Observations such as concentration measurements and flows, specified
  via
  [`sewer_data()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_data.md).
  The data can also be directly supplied to the relevant model
  components, in which case `data` can be left empty.

- assumptions:

  Assumptions about infections, shedding, or the sewage system,
  specified via
  [`sewer_assumptions()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_assumptions.md).
  The assumptions can also be directly supplied to the relevant model
  components, in which case `assumptions` can be left empty.

- measurements:

  The `measurements` module, see
  [`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md).

- sampling:

  The `sampling` module, see
  [`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md).

- sewage:

  The `sewage` module, see
  [`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md).

- shedding:

  The `shedding` module, see
  [`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md).

- infections:

  The `infections` module, see
  [`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md).

- fit_opts:

  Settings for model fitting, see
  [`set_fit_opts()`](https://adrian-lison.github.io/EpiSewer/reference/set_fit_opts.md).

- results_opts:

  Settings for results to be returned, see
  [`set_results_opts()`](https://adrian-lison.github.io/EpiSewer/reference/set_results_opts.md).

- run_fit:

  If `TRUE` (the default), the model is fitted immediately. If `FALSE`,
  the EpiSewerJob object is returned without fitting the model (it can
  be fitted later or even on a different machine).

## Value

If the model fitting was successful, a `list` with the following
elements is returned:

- `job`: the `EpiSewerJob` that was run

- `summary`: a summary of parameter estimates of interest

- `fitted`: the fitted model (if `results_opts(fitted = TRUE)`)

If the model fitting fails, a `list` with the following elements is
returned:

- `error`: information about the errors and warnings that were thrown

- `sampler_output`: potential outputs printed by the sampler
