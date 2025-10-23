# Estimate measurement noise

This option estimates the unexplained variation in wastewater
measurements using a constant coefficient of variation model.

If multiple measurements (replicates) per sample are provided,
`EpiSewer` can also explicitly model variation before the replication
stage.

For a non-constant coefficient of variation model, see
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md).

## Usage

``` r
noise_estimate(
  replicates = FALSE,
  cv_prior_mu = 0,
  cv_prior_sigma = 1,
  pre_replicate_cv_prior_mu = 0,
  pre_replicate_cv_prior_sigma = 1,
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

When `replicates=TRUE`, two coefficients of variation are estimated:

- the CV before the replication stage (see `pre_replicate_cv_prior_mu`)

- the CV after the replication stage (see `cv_prior_mu`)

The meaning of these CV estimates depends on the type of replicates. If
the replicates are biological replicates (i.e. independently processed),
then `cv` estimates the noise in the preprocessing and PCR, and
`pre_replicate_cv` estimates the noise from anything before
preprocessing (e.g. sampling noise and all other unexplained variation).
In contrast, if the replicates are technical replicates (i.e. several
PCR runs of the same preprocessed sample), then `cv` estimates only the
PCR noise, and `pre_replicate_cv` estimates all other noise (including
preprocessing noise.)

The priors of this component have the following functional form:

- coefficient of variation of concentration measurements (`cv`):
  `Truncated normal`

- coefficient of variation of concentration before the replication stage
  (`pre_replicate_cv`): `Truncated normal`

## See also

Other noise models:
[`noise_estimate_constant_var()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_constant_var.md),
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
