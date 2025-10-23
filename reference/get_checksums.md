# Get checksums that uniquely identify an EpiSewer job.

Get checksums that uniquely identify an EpiSewer job.

## Usage

``` r
get_checksums(job, stanmodel_instance = NULL)
```

## Arguments

- job:

  An EpiSewer job object.

- stanmodel_instance:

  A stanmodel instance.

## Value

A `list` with checksums identifying the model, input, inits, fit_opts
and results_opts.
