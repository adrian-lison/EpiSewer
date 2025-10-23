# Estimate measurement noise with constant variance

This option estimates the unexplained variation in wastewater
measurements using a constant variance model. This is usually a
misspecification and is only supported for comparison purposes.

For a constant coefficient of variation model, see
[`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md),
and for a non-constant coefficient of variation model, see
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md).

If multiple measurements (replicates) per sample are provided,
`EpiSewer` can also explicitly model variation before the replication
stage.

## Usage

``` r
noise_estimate_constant_var(
  replicates = FALSE,
  cv_prior_mu = 0,
  cv_prior_sigma = 1,
  pre_replicate_cv_prior_mu = 0,
  pre_replicate_cv_prior_sigma = 1,
  warn = TRUE,
  modeldata = modeldata_init()
)
```

## Arguments

- replicates:

  Should replicates be used to explicitly model variation before the
  replication stage?

- cv_prior_mu:

  Prior (mean) on the coefficient of variation of concentration
  measurements. Note that when `replicates=TRUE`, this is only the CV
  after the replication stage (see details for more explanation).

- cv_prior_sigma:

  Prior (standard deviation) on the coefficient of variation of
  concentration measurements.

- pre_replicate_cv_prior_mu:

  Prior (mean) on the coefficient of variation of concentrations
  *before* the replication stage.

- pre_replicate_cv_prior_sigma:

  Prior (standard deviation) on the coefficient of variation of
  concentrations *before* the replication stage.

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

Note that although this model keeps the variance constant, the prior for
the measurement noise is still in terms of the (average) coefficient of
variation (CV). This makes prior specification easier since the CV is
unitless.

The priors of this component have the following functional form:

- coefficient of variation of concentration measurements:
  `Truncated normal`

- coefficient of variation of concentration before the replication
  stage: `Truncated normal`

## See also

Other noise models:
[`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md),
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
