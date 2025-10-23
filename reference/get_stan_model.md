# Load stan model

Loads a specific stan model for model fitting in `EpiSewer`.

## Usage

``` r
get_stan_model(
  model_metainfo = NULL,
  model_filename = NULL,
  model_folder = "stan",
  profile = TRUE,
  threads = FALSE,
  force_recompile = FALSE,
  package = "EpiSewer"
)
```

## Arguments

- model_metainfo:

  A `list` containing meta information about the specified model to be
  fitted. The required stan model is automatically inferred from the
  `model_metainfo`.

- model_filename:

  File name of a specific stan model to load. This is an alternative to
  supplying `model_metainfo`.

- model_folder:

  Path to the folder containing the stan models for `EpiSewer`.

- profile:

  Should profiling be run during model fitting? Default is `TRUE`. If
  `FALSE`, will remove all profiling statements from the model before
  fitting. This can decrease runtime in some cases.

- threads:

  Should multihreading be enabled? Default is `FALSE`, as `EpiSewer`
  currently does not support within-chain parallelism.

- force_recompile:

  If `FALSE` (default), the model is only recompiled if changes to the
  model code are detected. However, as the change detection is not fully
  reliable, it is sometimes necessary to force recompilation after
  having made changes to the stan code.

- package:

  Name of the package in which to search for stan files (defaults to
  EpiSewer). If NULL, will search in the normal working directory.

## Value

A `list` containing the model definition and a link to the compiled stan
model.
