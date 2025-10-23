# Provide initialization value for a parameter based on the supplied prior with location and scale

Initialization using the prior is often better than initializing with
zero (and if the parameter is strictly positive, zero is not possible at
all). This function provides as init value the location of the prior
plus 1/4 of the scale. This ensure a positive init even if the mean is
zero (useful for truncated normal priors for example.)

## Usage

``` r
init_from_location_scale_prior(prior)
```

## Arguments

- prior:

  Prior for parameter as provided by
  [`set_prior()`](https://adrian-lison.github.io/EpiSewer/reference/set_prior.md).
  Should be a location and scale prior (first element location, second
  element scale).

## Value

Location of the prior plus 1/4 of the scale.

## Details

If the provided prior has zero variance, it is assumed that the
parameter will not be sampled and an empty init is returned.
