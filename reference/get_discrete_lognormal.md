# Get PMF of a discretized lognormal distribution

This function accepts both log-scale and unit-scale parameters to
specify a discretized lognormal distribution.

## Usage

``` r
get_discrete_lognormal(
  meanlog,
  sdlog,
  unit_mean = NULL,
  unit_sd = NULL,
  unit_cv = NULL,
  maxX = NULL,
  include_zero = TRUE,
  print_params = FALSE
)
```

## Arguments

- meanlog:

  Log scale mean (location of lognormal).

- sdlog:

  Log scale standard deviation (scale of lognormal).

- unit_mean:

  Alternative parameterization: unit/natural scale mean.

- unit_sd:

  Alternative parameterization: unit/natural scale sd.

- unit_cv:

  Alternative parameterization: unit/natural scale coefficient of
  variation.

- maxX:

  Right truncation point. All probability mass beyond `maxX` will be
  assigned to `maxX`. If `NULL` (default), this is automatically chosen
  such that the last bin has less than 0.5% of the probability mass.

- include_zero:

  Should the distribution explicitly cover X=0, or should X=1 include
  the probability mass for X=0 too?

- print_params:

  Should the log-level parameters be printed?

## Value

A numeric vector representing the PMF of the discretized lognormal.
