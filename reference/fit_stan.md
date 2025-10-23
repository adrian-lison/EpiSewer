# Fit model via stan

Helper function to fit an EpiSewer model using one of stan's algorithms.

## Usage

``` r
fit_stan(stanmodel_instance, arguments, fit_method, silent = FALSE)
```

## Arguments

- stanmodel_instance:

  A CmdStanModel object of the EpiSewer model constructed using
  [`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).

- arguments:

  Arguments to be passed to stan, including the data and the sampler
  options.

- fit_method:

  The algorithm used for sampling, currently supported are "mcmc" and
  "pathfinder".

- silent:

  Should the sampling be completely silent or should status updates and
  warnings be shown?

## Value

Either a fitted stan model with draws, or a list with error messages and
sampler outputs.
