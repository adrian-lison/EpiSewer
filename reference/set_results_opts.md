# Configure results returned after model fitting

Configure results returned after model fitting

## Usage

``` r
set_results_opts(
  fitted = TRUE,
  summary_intervals = c(0.5, 0.95),
  samples_ndraws = 50
)
```

## Arguments

- fitted:

  If `TRUE` (default), the fitted model object is also returned, not
  only summaries of the model fitting. Note that this makes the results
  object substantially larger.

- summary_intervals:

  Which credible intervals (CrIs) should be used to summarize the
  posterior distributions? Default is `c(0.5, 0.95)`, i.e. the 50% and
  95% CrI are returned

- samples_ndraws:

  Number of exemplary posterior samples to return. Note that the
  summaries always use all available samples.

## Value

A `list` with settings for the results.

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

ww_results_opts <- set_results_opts(
  fitted = TRUE, # return full model fit
  summary_intervals = c(0.5, 0.8, 0.95),
  samples_ndraws = 100
)

if (FALSE) ww_result <- EpiSewer(
  data = ww_data,
  assumptions = ww_assumptions,
  fit_opts = ww_fit_opts,
  results_opts = ww_results_opts
) # \dontrun{}
```
