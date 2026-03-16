# Estimate measurement noise for digital PCR data

This option models concentration measurements using an error model
specialized for digital PCR (e.g. digital droplet PCR). This is a fast
approximation assuming a fixed number of valid total partitions in the
dPCR assay and a known conversion factor. For a model that jointly
estimates uncertain assay parameters, see
[`noise_estimate_dPCR_params()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR_params.md).

This model predicts more relative variation at smaller concentrations,
which often leads to a better model fit, even for measurements from
other quantification methods such as qPCR.

If multiple measurements (replicates) per sample are provided,
`EpiSewer` can also explicitly model variation before the replication
stage.

## Usage

``` r
noise_estimate_dPCR(
  replicates = FALSE,
  cv_prior_mu = 0,
  cv_prior_sigma = 1,
  total_partitions_observe = FALSE,
  max_partitions = 30000,
  partition_loss_mean = 0.1,
  volume_scaled = 1e-05,
  pre_replicate_cv_prior_mu = 0,
  pre_replicate_cv_prior_sigma = 1,
  prePCR_noise_type = "log-normal",
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

- max_partitions:

  The maximum total number of dPCR partitions. This is usually defined
  by the manufacturer of the dPCR system/chip used, which supports a
  certain maximum number of partitions.

- partition_loss_mean:

  The mean relative partition loss. A certain proportion of partitions
  in a dPCR run is typically invalid and discarded from the
  concentration estimate.

- volume_scaled:

  The conversion factor (partition volume multiplied with the scaling of
  concentration in the assay) for the dPCR reaction. See details for
  further explanation.

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

The conversion factor (see `volume_scaled`) is the partition volume v
multiplied with a scaling factor s. The scaling factor accounts for
concentration differences between the sample and the reaction mix, for
example due to extraction or adding of reagents. For example, if the
partition volume is 4.5e-7 mL and the scaling factor is 100:3 (i.e. 100
gc/mL in the original sample correspond to 3 gc/mL in the PCR reaction),
then the overall conversion factor is 4.5e-7 \* 100 / 3 = 1.5e-5.

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

- coefficient of variation of concentration before the replication stage
  (`pre_replicate_cv`): `Truncated normal`

## See also

[LOD_estimate_dPCR](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
for a limit of detection model specialised for dPCR.

[noise_estimate_dPCR_params](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR_params.md)
for jointly estimating uncertain parameters of the assay.

Other noise models:
[`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md),
[`noise_estimate_constant_var()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_constant_var.md),
[`noise_estimate_dPCR_params()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR_params.md)
