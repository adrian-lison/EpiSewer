# Use the stan MCMC sampler

This option uses stan's NUTS sampler via
[cmdstanr::cmdstanr](https://mc-stan.org/cmdstanr/reference/cmdstanr-package.html)
for Markov Chain Monte Carlo (MCMC) sampling of the `EpiSewer` model.

## Usage

``` r
sampler_stan_mcmc(
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.99,
  max_treedepth = 15,
  step_size = 0.01,
  parallel_chains = NULL,
  threads_per_chain = 1,
  seed = 0,
  refresh = 200,
  show_messages = TRUE,
  show_exceptions = FALSE,
  init_pathfinder = TRUE,
  init_pathfinder_max_lbfgs_iters = NULL,
  ...
)
```

## Arguments

- chains:

  (positive integer) The number of Markov chains to run. The default is
  4.

- iter_warmup:

  (positive integer) The number of warmup iterations to run per chain.
  Note: in the CmdStan User's Guide this is referred to as `num_warmup`.

- iter_sampling:

  (positive integer) The number of post-warmup iterations to run per
  chain. Note: in the CmdStan User's Guide this is referred to as
  `num_samples`.

- adapt_delta:

  (real in `(0,1)`) The adaptation target acceptance statistic.

- max_treedepth:

  (positive integer) The maximum allowed tree depth for the NUTS engine.
  See the *Tree Depth* section of the CmdStan User's Guide for more
  details.

- step_size:

  (positive real) The *initial* step size for the discrete approximation
  to continuous Hamiltonian dynamics. This is further tuned during
  warmup.

- parallel_chains:

  (positive integer) The *maximum* number of MCMC chains to run in
  parallel. If `parallel_chains` is not specified then the default is to
  look for the option `"mc.cores"`, which can be set for an entire R
  session by `options(mc.cores=value)`. If the `"mc.cores"` option has
  not been set then the default is `1`.

- threads_per_chain:

  (positive integer) If the model was
  [compiled](https://mc-stan.org/cmdstanr/reference/model-method-compile.html)
  with threading support, the number of threads to use in parallelized
  sections *within* an MCMC chain (e.g., when using the Stan functions
  `reduce_sum()` or `map_rect()`). This is in contrast with
  `parallel_chains`, which specifies the number of chains to run in
  parallel. The actual number of CPU cores used is
  `parallel_chains*threads_per_chain`. For an example of using threading
  see the Stan case study [Reduce Sum: A Minimal
  Example](https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html).

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

- init_pathfinder:

  Should stan's pathfinder algorihtm be used to initialize the MCMC
  chains?

- init_pathfinder_max_lbfgs_iters:

  The maximum number of iterations for LBFGS during pathfinder
  initialization.

- ...:

  Further arguments to pass to
  [cmdstanr::cmdstanr](https://mc-stan.org/cmdstanr/reference/cmdstanr-package.html).

## Value

A `list` with settings for the MCMC sampler.
