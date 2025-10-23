# Configure the model fitting

Sets model fitting options in `EpiSewer`.

## Usage

``` r
set_fit_opts(sampler = sampler_stan_mcmc(), model = model_stan_opts())
```

## Arguments

- sampler:

  Which type of sampler should be used for model fitting? Currently
  supported are

  - [`sampler_stan_mcmc()`](https://adrian-lison.github.io/EpiSewer/reference/sampler_stan_mcmc.md):
    default, recommended for inference

  - [`sampler_stan_pathfinder()`](https://adrian-lison.github.io/EpiSewer/reference/sampler_stan_pathfinder.md):
    potentially inaccurate, recommended for preview only

- model:

  Details about the model file to be used, see
  [`model_stan_opts()`](https://adrian-lison.github.io/EpiSewer/reference/model_stan_opts.md).
  This also makes it possible for users to supply customized models.

## Value

A `list` with the model fitting options.

## Examples

``` r
ww_data <- ww_data_SARS_CoV_2_Zurich
ww_assumptions <- ww_assumptions_SARS_CoV_2_Zurich

ww_fit_opts <- set_fit_opts(
  sampler = sampler_stan_mcmc(
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 400,
    iter_sampling = 400,
    seed = 42 # ensures reproducibility
  )
)

if (FALSE) ww_result <- EpiSewer(
  data = ww_data,
  assumptions = ww_assumptions,
  fit_opts = ww_fit_opts
) # \dontrun{}
```
