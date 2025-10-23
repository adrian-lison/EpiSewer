# Run stan's pathfinder algorithm to obtain improved inits.

Helper function that runs the pathfinder algorithm with default inits
and returns the result to be used as inits in mcmc sampling.

## Usage

``` r
get_pathfinder_inits(stanmodel_instance, job, max_iters = NULL)
```

## Arguments

- stanmodel_instance:

  A CmdStanModel object of the EpiSewer model constructed using
  [`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).

- job:

  An EpiSewer job object with data and default inits.

- max_iters:

  The maximum number of iterations for LBFGS during pathfinder
  initialization.

## Value

Fitted model, a CmdStanPathfinder object. Can be used as arguments for
the inits of the `CmdStanModel$sample()` method.
