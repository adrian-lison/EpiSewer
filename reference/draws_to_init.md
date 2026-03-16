# Convert draws from a previous stan model fit to inits for MCMC sampling

Helper function that converts the draws obtained (e.g. from pathfinder)
to a list of inits that can be used for MCMC sampling.

## Usage

``` r
draws_to_init(fit, model, num_chains = 4)
```

## Arguments

- fit:

  A fitted model object (e.g. obtained from pathfinder).

- model:

  The model specification, used to identify scalar parameters.

- num_chains:

  The number of chains for which to generate inits.

## Value

A list of inits for MCMC sampling, with one entry per chain.
