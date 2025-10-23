# Model the infection process

This module function is used to specify the components of the
`infections` module in `EpiSewer`.

Each component can be specified using one or several helper functions
(see available options below). See the documentation of the individual
helper functions to adjust model priors and further settings.

## Usage

``` r
model_infections(
  generation_dist = generation_dist_assume(),
  R = R_estimate_gp(),
  seeding = seeding_estimate_rw(),
  infection_noise = infection_noise_estimate()
)
```

## Arguments

- generation_dist:

  Generation time distribution. The intrinsic distribution of the time
  between infection of a primary case and infection of its secondary
  cases. Modeling options:

  - [`generation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/generation_dist_assume.md)

- R:

  Effective reproduction number over time. This is the main parameter of
  interest estimated by `EpiSewer`. `R` is smoothed using a time series
  smoothing prior (Gaussian process by default). A variety of smoothing
  priors is supported. Modeling options:

  - [`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md)

  - [`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md)

  - [`R_estimate_smooth_derivative()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_smooth_derivative.md)

  - [`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md)

  - [`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md)

  - [`R_estimate_approx()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_approx.md)

  - [`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md)

  - [`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)

- seeding:

  Seeding of initial infections. The renewal model used by `EpiSewer`
  requires a seeding phase of the length of the maximum generation time.
  For these initial infections, a simple seeding model instead of the
  renewal model must be used. Modeling options:

  - [`seeding_estimate_growth()`](https://adrian-lison.github.io/EpiSewer/reference/seeding_estimate_growth.md)

  - [`seeding_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/seeding_estimate_rw.md)

  - [`seeding_estimate_constant()`](https://adrian-lison.github.io/EpiSewer/reference/seeding_estimate_constant.md)

- infection_noise:

  Noise in the infection process. `EpiSewer` implements a stochastic
  infection model, i.e. allows for variation in the number of new
  infections generated at each time step. This accounts for stochastic
  uncertainty in the infection process and often speeds up model
  fitting. Modeling options:

  - [`infection_noise_none()`](https://adrian-lison.github.io/EpiSewer/reference/infection_noise_none.md)

  - [`infection_noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/infection_noise_estimate.md)

## Value

A `modeldata` object containing the data and specifications of the
`infections` module.

## See also

Other module functions:
[`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md),
[`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md),
[`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md),
[`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md),
[`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md)
