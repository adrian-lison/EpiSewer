# Define a normal prior in modeldata

Define a normal prior in modeldata

## Usage

``` r
set_prior_normal(
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

  Mean.

- sigma:

  Standard deviation.

- two_sigma:

  Two times the standard deviation. Useful for defining priors via the
  two-sigma rule-of-thumb (approximately 95% of probability mass).

- q5:

  Lower quantile (5%).

- q95:

  Upper quantile (95%).

## Value

Prior specification for modeldata.
