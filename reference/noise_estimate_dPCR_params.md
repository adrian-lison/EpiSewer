# Estimate measurement noise for digital PCR data with priors for assay parameters

This option models concentration measurements using an error model
specialized for digital PCR (e.g. digital droplet PCR). This option also
jointly estimates the uncertain parameters of the dPCR assay. For a
faster approximation assuming fixed assay parameters, see
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md).

This model predicts more relative variation at smaller concentrations,
which often leads to a better model fit, even for measurements from
other quantification methods such as qPCR.

If multiple measurements (replicates) per sample are provided,
`EpiSewer` can also explicitly model variation before the replication
stage.

## Usage

``` r
noise_estimate_dPCR_params(
  replicates = FALSE,
  cv_prior_mu = 0,
  cv_prior_sigma = 1,
  total_partitions_observe = FALSE,
  max_partitions_prior_lower = 20000,
  max_partitions_prior_upper = 40000,
  partition_loss_mean_prior_lower = 0.05,
  partition_loss_mean_prior_upper = 0.3,
  partition_loss_variation_prior_lower = 0.5,
  partition_loss_variation_prior_upper = 2,
  partition_loss_max = 0.5,
  volume_scaled_prior_lower = 1e-06,
  volume_scaled_prior_upper = 1e-04,
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

  Prior (mean) on the coefficient of variation of analyzed
  concentrations (pre-PCR noise). Note that this does *not* include the
  technical noise of the digital PCR. This is because the dPCR noise is
  explicitly modeled (using the total number of partitions and
  conversion factor). Moreover, when `replicates=TRUE`, noise before the
  replication stage is separately modeled, so `cv_prior_mu` only
  describes any remaining pre-PCR noise (see details for more
  explanation).

- cv_prior_sigma:

  Prior (standard deviation) on the coefficient of variation of analyzed
  concentrations.

- total_partitions_observe:

  If TRUE, the total number of partitions is taken from the supplied
  measurements `data.frame`. This requires that the argument
  `total_partitions_col` is specified in
  [`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md).

- max_partitions_prior_lower:

  Prior (lower bound) for the maximum total number of dPCR partitions.
  This is usually defined by the manufacturer of the dPCR system/chip
  used, which supports a certain maximum number of partitions. If you
  know the exact dPCR system and its maximum partition number, you can
  set both `max_partitions_prior_lower` and `max_partitions_prior_upper`
  to this value. Otherwise, this prior can be used to set a broad lower
  and upper bound for the maximum number of partitions, to reflect a
  range of popular dPCR systems/chips.

- max_partitions_prior_upper:

  Prior (upper bound) for the maximum total number of dPCR partitions
  (see `max_partitions_prior_lower` for details.)

- partition_loss_mean_prior_lower:

  Prior (5% quantile) for the mean relative partition loss. A certain
  proportion of partitions in a dPCR run is typically invalid and
  discarded from the concentration estimate. This prior can be used to
  set a lower and upper bound for the mean proportion of partitions
  lost. Note that for proportions close to 0 or to `partition_loss_max`
  (see below), the resulting mean partition loss can differ from what is
  specified here, because we internally translate this prior to the
  logit scale.

- partition_loss_mean_prior_upper:

  Prior (95% quantile) for the mean relative partition loss (see
  `partition_loss_mean_prior_lower` for details). In well-functioning
  dPCR assays, the mean proportion of lost partitions should not be very
  high (definitely below 50%).

- partition_loss_variation_prior_lower:

  Prior (5% quantile) for the variation in the number of invalid
  partitions across dPCR runs. The proportion of partitions lost
  typically varies between dPCR runs. We thus model this proportion as
  logit-normal distributed, with mean (approximately) defined by
  `partition_loss_mean_prior` and logit-level standard deviation
  `sigma`. You can use `partition_loss_variation_prior_lower` and
  `partition_loss_variation_prior_upper` to set a lower and upper bound
  for `sigma`.

- partition_loss_variation_prior_upper:

  Prior (95% quantile) for the variation in the number of invalid
  partitions across dPCR runs (see
  `partition_loss_variation_prior_lower` for details).

- partition_loss_max:

  The maximum proportion of partitions that can be lost in a valid dPCR
  run. During quality control, runs where the proportion of invalid
  partitions is above some threshold (e.g. 50%) are often discarded.
  This parameter can be used to represent such a QC threshold.

- volume_scaled_prior_lower:

  Prior (lower bound) on the conversion factor (partition volume scaled
  by the dilution of wastewater in the assay) for the dPCR reaction. See
  details for further explanation.

- volume_scaled_prior_upper:

  Prior (upper bound) on the conversion factor (partition volume scaled
  by the dilution of wastewater in the assay) for the dPCR reaction. If
  this is identical to `volume_scaled_prior_lower`, the conversion
  factor will be fixed and not estimated.

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

The conversion factor (see `volume_scaled_prior_lower`,
`volume_scaled_prior_upper`) is the partition volume v multiplied with a
scaling factor s. The scaling factor accounts for concentration
differences between the sample and the reaction mix, for example due to
extraction or adding of reagents. For example, if the partition volume
is 4.5e-7 mL and the scaling factor is 100:3 (i.e. 100 gc/mL in the
original sample correspond to 3 gc/mL in the PCR reaction), then the
overall conversion factor is 4.5e-7 \* 100 / 3 = 1.5e-5.

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

- maximum number of total partitions in dPCR: `Uniform`

- mean proportion of lost partitions in dPCR: `Normal (logit-level)`

- variation of proportion of lost partitions:
  `Truncated normal (logit-level)`

- conversion factor for dPCR: `Uniform`

- coefficient of variation of concentration before the replication stage
  (`pre_replicate_cv`): `Truncated normal`

## See also

[LOD_estimate_dPCR](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
for a limit of detection model specialised for dPCR.

[noise_estimate_dPCR](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
for a faster dPCR noise model with fixed assay parameters.

Other noise models:
[`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md),
[`noise_estimate_constant_var()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_constant_var.md),
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
