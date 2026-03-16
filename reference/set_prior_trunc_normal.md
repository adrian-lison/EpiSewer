# Define a truncated normal prior in modeldata

This function defines a truncated normal prior (truncation at zero) for
a parameter in modeldata.

## Usage

``` r
set_prior_trunc_normal(
  param,
  mu = NULL,
  sigma = NULL,
  two_sigma = NULL,
  q5 = NULL,
  q95 = NULL
)
```

## Arguments

- param:

  Name of the parameter for which the prior is defined.

- mu:

  Mean (not accounting for truncation).

- sigma:

  Standard deviation (not accounting for truncation).

- two_sigma:

  Two times the standard deviation (not accounting for truncation).
  Useful for defining priors via the two-sigma rule-of-thumb
  (approximately 95% of probability mass).

- q5:

  Lower quantile (5%).

- q95:

  Upper quantile (95%).

## Value

Prior specification for modeldata.
