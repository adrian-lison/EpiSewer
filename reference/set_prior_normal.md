# Define a normal prior in modeldata

Define a normal prior in modeldata

## Usage

``` r
set_prior_normal(param, mu, sigma = NULL, two_sigma = NULL)
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

## Value

Prior specification for modeldata.
