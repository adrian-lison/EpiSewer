# Find mu and sigma of a truncated normal distribution for given quantiles

This function estimates the parameters `mu` and `sigma` of a normal
distribution truncated at zero, such that the 5% and 95% quantiles of
the distribution match the specified target quantiles.

## Usage

``` r
find_mu_sigma_truncnorm(q5_target, q95_target)
```

## Arguments

- q5_target:

  The target quantile at 5% (lower bound).

- q95_target:

  The target quantile at 95% (upper bound).

## Value

A list containing the estimated `mu` and `sigma` for the truncated
normal distribution that matches the specified quantiles.
