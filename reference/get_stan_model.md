# Load stan model

Loads a specific stan model for model fitting in `EpiSewer`.

## Usage

``` r
get_stan_model(
  model_filename = NULL,
  model_folder = "stan",
  profile = TRUE,
  threads = FALSE,
  force_recompile = FALSE,
  package = "EpiSewer",
  use_docker = FALSE
)
```

## Arguments

- model_filename:

  File name of a specific stan model to load. If NULL, the default model
  will be loaded.

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

- use_docker:

  Use the model compiled in a docker container? Note that if TRUE, all
  other arguments are ignored and the precompiled model in the container
  will be used. Default is FALSE.

## Value

A `list` containing the model definition and a link to the compiled stan
model.
