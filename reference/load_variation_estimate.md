# Estimate individual-level load variation

This option accounts for variation in the total shedding load per case
by modeling individual shedding loads as Gamma distributed with mean
equal to the average `load_per_case` and a coefficient of variation to
be estimated.

## Usage

``` r
load_variation_estimate(
  cv_prior_mu = 1,
  cv_prior_sigma = 0,
  modeldata = modeldata_init()
)
```

## Arguments

- cv_prior_mu:

  Mean of the truncated normal prior for the coefficient of
  individual-level variation. Default is a CV of 1, which means that
  approximately 60% of individuals shed less than the mean, and
  approximately 10% of individuals shed less than 10% of the mean.

- cv_prior_sigma:

  Standard deviation of the truncated normal prior for the coefficient
  of individual-level variation. If this is set to zero, the coefficient
  of variation is fixed to the mean.

- modeldata:

  A `modeldata` object to which the above model specifications should be
  added. Default is an empty model given by
  [`modeldata_init()`](https://adrian-lison.github.io/EpiSewer/reference/modeldata_init.md).
  Can also be an already partly specified model returned by other
  `EpiSewer` modeling functions.

## Value

A `modeldata` object containing data and specifications of the model to
be fitted. Can be passed on to other `EpiSewer` modeling functions to
add further data and model specifications.

The `modeldata` object also includes information about parameter
initialization (`.init`), meta data (`.metainfo`), and checks to be
performed before model fitting (`.checks`).

## Details

Note that measurement noise and individual-level load variation might
not be jointly identifiable. This is why the coefficient of variation of
the individual-level load is fixed by default.

Also note that the accuracy of the variation model depends on the
estimated number of cases to be on the right scale (which in turn
depends on the assumed `load_per_case` to be roughly correct). This is
because in the variation model, the population-level coefficient of
variation (CV) of shedding loads is proportional to the individual-level
CV divided by the square root of the number of cases. If for example the
assumed `load_per_case` is too small, the number of cases will be
overestimated, which has effects in two directions:

- the individual-level CV will be overestimated (especially if the prior
  on the individual-level CV is weak)

- the population-level CV will be underestimated (especially if the
  prior on the individual-level CV is strong)

## See also

Other load variation models:
[`load_variation_none()`](https://adrian-lison.github.io/EpiSewer/reference/load_variation_none.md)
