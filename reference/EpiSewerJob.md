# Constructor for EpiSewerJob objects

Constructor for EpiSewerJob objects

## Usage

``` r
EpiSewerJob(
  job_name,
  modeldata,
  fit_opts,
  results_opts,
  jobarray_size = 1,
  overwrite = TRUE,
  results_exclude = c()
)
```

## Arguments

- job_name:

  A name for the job, used for labeling outputs and results.

- modeldata:

  A list with all model specifications.

- fit_opts:

  Settings for model fitting, see
  [`set_fit_opts()`](https://adrian-lison.github.io/EpiSewer/reference/set_fit_opts.md).

- results_opts:

  Settings for results to be returned, see
  [`set_results_opts()`](https://adrian-lison.github.io/EpiSewer/reference/set_results_opts.md).

- jobarray_size:

  The size of the job array to run. This is relevant when running the
  model on a cluster with job arrays. The default is 1, i.e. no job
  array.

- overwrite:

  If `TRUE` (the default), existing outputs with the same job name will
  be overwritten.

- results_exclude:

  A character vector with names of results to exclude when saving
  results. This can be used to exclude large objects such as the fitted
  model from being saved in the results.

## Value

An `EpiSewerJob` object with all specifications for model fitting.
