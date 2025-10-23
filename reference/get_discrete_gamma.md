# Get PMF of a discretized Gamma distribution

This function accepts different parameterizations to specify a
discretized Gamma distribution.

## Usage

``` r
get_discrete_gamma(
  gamma_shape,
  gamma_rate,
  gamma_scale,
  gamma_mean,
  gamma_sd,
  gamma_cv,
  maxX = NULL,
  include_zero = TRUE,
  print_params = FALSE
)
```

## Arguments

- gamma_shape:

  Shape parameter of the Gamma distribution.

- gamma_rate:

  Rate parameter of the Gamma distribution.

- gamma_scale:

  Scale parameter of the Gamma distribution. Can be specified instead of
  the rate. Only has an effect if the rate is not specified.

- gamma_mean:

  Alternative parameterization: Mean of the Gamma distribution.

- gamma_sd:

  Alternative parameterization: Standard deviation of the Gamma
  distribution.

- gamma_cv:

  Alternative parameterization: Coefficient of variation of the Gamma
  distribution.

- maxX:

  Right truncation point. All probability mass beyond `maxX` will be
  assigned to `maxX`. If `NULL` (default), this is automatically chosen
  such that the last bin has less than 0.5% of the probability mass.

- include_zero:

  Should the distribution explicitly cover X=0, or should X=1 include
  the probability mass for X=0 too?

- print_params:

  Should the shape and rate parameters be printed?

## Value

A numeric vector representing the PMF of the discretized Gamma
distribution.
