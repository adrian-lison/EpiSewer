# Estimate sample effects using a design matrix

This option uses a log-linear regression model to estimate effects of
sample covariates on the concentration. Concentrations can be influenced
by sampling-related external factors, for example the time between
sampling and shipping to the lab (age-of-sample effect), or different
sampling or storage methods.

## Usage

``` r
sample_effects_estimate_matrix(
  design_matrix,
  effect_prior_mu = 0,
  effect_prior_sigma = 1,
  modeldata = modeldata_init()
)
```

## Arguments

- design_matrix:

  A design matrix with different covariates potentially influencing
  sample concentration. The matrix must have one row for each modeled
  day. See
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  for construction of design matrices.

- effect_prior_mu:

  Prior (mean) on the regression coefficients.

- effect_prior_sigma:

  Prior (standard deviation) on the regression coefficients.

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

`EpiSewer` will fit a fixed-effects log-linear model, random effects are
currently not supported.

The priors of this component have the following functional form:

- regression coefficients: `Normal`

## See also

Other sample effect models:
[`sample_effects_estimate_weekday()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_weekday.md),
[`sample_effects_none()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_none.md)
