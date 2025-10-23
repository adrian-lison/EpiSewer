# Estimate measurement noise for digital PCR data

This option estimates the unexplained variation in wastewater
measurements using a coefficient of variation model specialized for
digital PCR (e.g. digital droplet PCR). Specifically, the coefficient of
variation is modeled as a function of the expected concentration
according to the statistical properties of dPCR.

This model predicts a higher coefficient of variation at smaller
concentrations, which often leads to a better model fit, even for
measurements from other quantification methods such as qPCR.

If multiple measurements (replicates) per sample are provided,
`EpiSewer` can also explicitly model variation before the replication
stage.

## Usage

``` r
noise_estimate_dPCR(
  replicates = FALSE,
  cv_prior_mu = 0,
  cv_prior_sigma = 1,
  total_partitions_prior_mu = 20000,
  total_partitions_prior_sigma = 5000,
  total_partitions_observe = FALSE,
  partition_variation_prior_mu = 0,
  partition_variation_prior_sigma = 0.05,
  volume_scaled_prior_mu = 1e-05,
  volume_scaled_prior_sigma = 4e-05,
  pre_replicate_cv_prior_mu = 0,
  pre_replicate_cv_prior_sigma = 1,
  prePCR_noise_type = "log-normal",
  use_taylor_approx = TRUE,
  modeldata = modeldata_init()
)
```

## Arguments

- replicates:

  Should replicates be used to explicitly model variation before the
  replication stage?

- cv_prior_mu:

  Prior (mean) on the coefficient of variation of concentration
  measurements. Note that in contrast to using
  [`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md),
  this does *not* include the technical noise of the digital PCR. This
  is because the dPCR noise is explicitly modeled (using the number of
  partitions and conversion factor). Moreover, when `replicates=TRUE`,
  this is only the CV after the replication stage (see details for more
  explanation).

- cv_prior_sigma:

  Prior (standard deviation) on the coefficient of variation of
  concentration measurements.

- total_partitions_prior_mu:

  Prior (mean) on the total number of partitions in the dPCR reaction.

- total_partitions_prior_sigma:

  Prior (standard deviation) on the total number of partitions in the
  dPCR reaction. If this is set to zero, the total number of partitions
  will be fixed to the prior mean and not estimated.

- total_partitions_observe:

  If TRUE, the total number of partitions is taken from the supplied
  measurements `data.frame`. This requires that the argument
  `total_partitions_col` is specified in
  [`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md).

- partition_variation_prior_mu:

  Prior (mean) on the coefficient of variation of the total number of
  partitions in the dPCR reaction. Usually, the maximum number of
  partitions possible for a given dPCR chip is not reached, i.e. a
  certain number of partitions is lost. This loss varies between PCR
  runs, and is modeled as log-normal distributed in EpiSewer.

- partition_variation_prior_sigma:

  Prior (standard deviation) on the coefficient of variation of the
  total number of partitions in the dPCR reaction. If this is set to
  zero, the partition variation will be fixed to the prior mean and not
  estimated.

- volume_scaled_prior_mu:

  Prior (mean) on the conversion factor (partition volume scaled by the
  dilution of wastewater in the assay) for the dPCR reaction. See
  details for further explanation.

- volume_scaled_prior_sigma:

  Prior (standard deviation) on the conversion factor (partition volume
  scaled by the dilution of wastewater in the assay) for the dPCR
  reaction. If this is set to zero, the conversion factor will be fixed
  to the prior mean and not estimated.

- pre_replicate_cv_prior_mu:

  Prior (mean) on the coefficient of variation of concentrations
  *before* the replication stage.

- pre_replicate_cv_prior_sigma:

  Prior (standard deviation) on the coefficient of variation of
  concentrations *before* the replication stage.

- prePCR_noise_type:

  The parametric distribution to assume for noise before the PCR assay.
  Currently supported are "log-normal" and "gamma". The choice of the
  parametric distribution typically makes no relevant difference for the
  noise model, but can make a relevant difference for the LOD model if
  [`LOD_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
  is used.

- use_taylor_approx:

  If TRUE (default), a Taylor expansion approximation is used to
  estimate the CV of measurements under pre-PCR noise. The approximation
  is very accurate, unless concentrations are extremely high (so high
  that the quality of the measurements from dPCR would anyway be
  questionable).

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

The conversion factor (see `volume_scaled_prior_mu`,
`volume_scaled_prior_sigma`) is the partition volume scaled by the
dilution of the wastewater in the assay. The dilution accounts for all
extraction and preparation steps. For example, if the partition volume
is 4.5e-7 mL and the dilution of the wastewater is 100:3 (i.e. 100 gc/mL
in the original wastewater sample correspond to 3 gc/mL in the PCR
reaction), then the overall conversion factor is 4.5e-7 \* 100 / 3 =
1.5e-5.

When `replicates=TRUE`, two coefficients of variation are estimated:

- the CV before the replication stage (see `pre_replicate_cv_prior_mu`)

- the CV after the replication stage (see `cv_prior_mu`)

  The meaning of these CV estimates depends on the type of replicates.
  If the replicates are biological replicates (i.e. independently
  processed), then `cv` estimates the noise in the preprocessing before
  the PCR, and `pre_replicate_cv` estimates the noise from anything
  before preprocessing (e.g. sampling noise and all other unexplained
  variation). In contrast, if the replicates are technical replicates
  (i.e. several PCR runs of the same preprocessed sample), then `cv`
  estimates only unexplained PCR noise (should be close to zero), and
  `pre_replicate_cv` estimates all other noise (including preprocessing
  noise.)

The priors of this component have the following functional form:

- coefficient of variation of concentration measurements (`cv`):
  `Truncated normal`

- mean number of partitions in dPCR: `Truncated normal`

- coefficient of variation of number of partitions in dPCR:
  `Truncated normal`

- conversion factor for dPCR: `Truncated normal`

- coefficient of variation of concentration before the replication stage
  (`pre_replicate_cv`): `Truncated normal`

## See also

[LOD_estimate_dPCR](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
for a limit of detection model specialised for dPCR.

Other noise models:
[`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md),
[`noise_estimate_constant_var()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_constant_var.md)
