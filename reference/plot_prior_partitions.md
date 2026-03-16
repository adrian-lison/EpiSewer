# Plot prior distributions for total number of partitions

Plot prior distributions for total number of partitions

## Usage

``` r
plot_prior_partitions(modeldata, n_draws = 1000, show_draws = 50, seed = 0)
```

## Arguments

- modeldata:

  A `modeldata` object as returned by
  [`noise_estimate_dPCR_params()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR_params.md).

- n_draws:

  Number of draws to simulate from the prior distributions.

- show_draws:

  Number of example draws to show in the plots.

- seed:

  Seed for random number generation to ensure reproducibility.

## Value

A grid of plots showing the distribution of partition loss, 95%
intervals for valid partitions, and distributions of mean and standard
deviation of valid partitions.
