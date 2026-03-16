# Fit the EpiSewer model

Helper function that fits the EpiSewer model using the specified sampler
and model instance, with the data and inits from the job object.

## Usage

``` r
fit_model(job, model, model_instance, run_silent = FALSE)
```

## Arguments

- job:

  An EpiSewer job object containing the data, inits, and fitting
  options.

- model:

  The model specification, used to identify parameters.

- model_instance:

  A compiled model object.

- run_silent:

  Should the fitting be run silently? Default is FALSE.

## Value

A fitted model object if successful, or a list with error messages and
sampler output if there was an error during fitting.
