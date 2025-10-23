# Get PMF of a discretized shifted Gamma distribution

This function specifies a discretized shifted Gamma distribution
(shifted such that the minimum is at 1) using a mean and standard
deviation parameterization.

The shift makes the distribution attractive for modeling generation time
distributions, which are often assumed to have zero probability mass for
a generation time of zero (as this is incompatible with the renewal
equation).

## Usage

``` r
get_discrete_gamma_shifted(gamma_mean, gamma_sd, maxX = NULL)
```

## Arguments

- gamma_mean:

  Mean of the shifted Gamma distribution.

- gamma_sd:

  Standard deviation of the shifted Gamma distribution.

- maxX:

  Right truncation point. All probability mass beyond `maxX` will be
  assigned to `maxX`. If `NULL` (default), this is automatically chosen
  such that the last bin has approximately less than 0.5% of the
  probability mass.

## Value

A numeric vector representing the PMF of the discretized shifted Gamma
distribution.

## Details

This code was adapted from EpiEstim, credit to Anne Cori
(a.cori@imperial.ac.uk).
