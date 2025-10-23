# Use stan's pathfinder variational inference algorithm

This option uses stan's pathfinder algorithm via
[cmdstanr::cmdstanr](https://mc-stan.org/cmdstanr/reference/cmdstanr-package.html)
for variational inference. Sampling is fast but yields potentially
inaccurate results with unreliable uncertainty information. Currently,
this sampler is recommended only for preview purposes, not for
real-world inference.

## Usage

``` r
sampler_stan_pathfinder(
  draws = 4000,
  num_paths = 4,
  single_path_draws = NULL,
  max_lbfgs_iters = NULL,
  num_elbo_draws = NULL,
  psis_resample = FALSE,
  seed = 0,
  refresh = 1000,
  show_messages = TRUE,
  show_exceptions = FALSE,
  ...
)
```

## Arguments

- draws:

  (positive integer) Number of draws to return after performing pareto
  smooted importance sampling (PSIS). This should be smaller than
  `single_path_draws * num_paths` (future versions of CmdStan will throw
  a warning).

- num_paths:

  (positive integer) Number of single pathfinders to run.

- single_path_draws:

  (positive integer) Number of draws a single pathfinder should return.
  The number of draws PSIS sampling samples from will be equal to
  `single_path_draws * num_paths`.

- max_lbfgs_iters:

  (positive integer) The maximum number of iterations for LBFGS.

- num_elbo_draws:

  (positive integer) Number of draws to make when calculating the ELBO
  of the approximation at each iteration of LBFGS.

- psis_resample:

  (logical) Whether to perform pareto smoothed importance sampling. If
  `TRUE`, the number of draws returned will be equal to `draws`. If
  `FALSE`, the number of draws returned will be equal to
  `single_path_draws * num_paths`.

- seed:

  (positive integer(s)) A seed for the (P)RNG to pass to CmdStan. In the
  case of multi-chain sampling the single `seed` will automatically be
  augmented by the the run (chain) ID so that each chain uses a
  different seed. The exception is the transformed data block, which
  defaults to using same seed for all chains so that the same data is
  generated for all chains if RNG functions are used. The only time
  `seed` should be specified as a vector (one element per chain) is if
  RNG functions are used in transformed data and the goal is to generate
  *different* data for each chain.

- refresh:

  (non-negative integer) The number of iterations between printed screen
  updates. If `refresh = 0`, only error messages will be printed.

- show_messages:

  (logical) When `TRUE` (the default), prints all output during the
  execution process, such as iteration numbers and elapsed times. If the
  output is silenced then the
  [`$output()`](https://mc-stan.org/cmdstanr/reference/fit-method-output.html)
  method of the resulting fit object can be used to display the silenced
  messages.

- show_exceptions:

  (logical) When `TRUE` (the default), prints all informational
  messages, for example rejection of the current proposal. Disable if
  you wish to silence these messages, but this is not usually
  recommended unless you are very confident that the model is correct up
  to numerical error. If the messages are silenced then the
  [`$output()`](https://mc-stan.org/cmdstanr/reference/fit-method-output.html)
  method of the resulting fit object can be used to display the silenced
  messages.

- ...:

  Further arguments to pass to
  [cmdstanr::cmdstanr](https://mc-stan.org/cmdstanr/reference/cmdstanr-package.html).

## Value

A `list` with settings for the pathfinder sampler.
