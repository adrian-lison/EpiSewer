# Specify details of the stan model

Specify details of the stan model

## Usage

``` r
model_stan_opts(
  model_filename = NULL,
  model_folder = "stan",
  profile = TRUE,
  threads = FALSE,
  force_recompile = FALSE,
  package = "EpiSewer"
)
```

## Arguments

- model_filename:

  Name of the stan model file. If NULL (default), the name of the model
  file is automatically inferred based on the model specification.

- model_folder:

  (Relative) path to the folder containing the stan model.

- profile:

  Should profiling be run during model fitting? Default is `TRUE`. If
  `FALSE`, will remove all profiling statements from the model before
  fitting. This can decrease runtime in some cases.

- threads:

  Should multihreading be enabled? Default is `FALSE`, as `EpiSewer`
  currently does not support within-chain parallelism.

- force_recompile:

  Should recompilation be forced before model fitting? If `FALSE`
  (default), the model is only recompiled when changes to the model code
  are detected. However, as the change detection is not fully reliable,
  it is sometimes necessary to force recompilation after having made
  changes to the stan code.

- package:

  Name of the package in which to search for stan files (defaults to
  EpiSewer). If NULL, will search in the normal working directory,
  allowing users to provide their own customized model files.

## Value

A list with details of the stan model
