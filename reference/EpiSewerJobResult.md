# Constructor for EpiSewerJobResult objects

Constructor for EpiSewerJobResult objects

## Usage

``` r
EpiSewerJobResult(fit_res, job, stan_model, checksums)
```

## Arguments

- fit_res:

  The result of model fitting, either a fitted model object or an error
  object if model fitting failed.

- job:

  The EpiSewerJob object that was run.

- stan_model:

  The stan model that was used for fitting.

- checksums:

  A list with checksums that uniquely identify the model, input, inits,
  fit_opts and results_opts of the job.

## Value

An EpiSewerJobResult object with the results of model fitting and
relevant information about the job and the model.
