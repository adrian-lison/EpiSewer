# Estimate measurement noise (internal helper function)

This option estimates the unexplained variation in wastewater
measurements. If multiple measurements (replicates) per sample are
provided, `EpiSewer` can also explicitly model variation before the
replication stage.

This helper function is called from other noise modeling functions.
[`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md)
is a constant coefficient of variation model,
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
is a noise model specialized for digital PCR (`cv_type = "dPCR"`), which
may however also work with other quantification methods such as qPCR,
and
[`noise_estimate_constant_var()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_constant_var.md)
is a constant variance model.

## Usage

``` r
noise_estimate_(
  replicates = FALSE,
  cv_prior_mu = 0,
  cv_prior_sigma = 1,
  cv_type = "constant",
  total_partitions_prior_mu = NULL,
  total_partitions_prior_sigma = NULL,
  total_partitions_observe = NULL,
  partition_variation_prior_mu = NULL,
  partition_variation_prior_sigma = NULL,
  volume_scaled_prior_mu = NULL,
  volume_scaled_prior_sigma = NULL,
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
  measurements. Note that when `replicates=TRUE`, this is only the CV
  after the replication stage (see details for more explanation).

- cv_prior_sigma:

  Prior (standard deviation) on the coefficient of variation of
  concentration measurements.

- cv_type:

  One out of "constant" (default), "constant_var", or "dPCR". If
  "constant", the coefficient of variation is estimated as a
  constant/single parameter for all observations. If "dPCR", the
  coefficient of variation is modeled as a function of the expected
  concentration according to the statistical properties of dPCR. In
  particular, this model predicts a higher coefficient of variation at
  smaller concentrations, which often leads to a better model fit. If
  "constant_var", not the coefficient of variation but the variance of
  measurements is modeled as constant. This is usually a
  misspecification and is only supported for comparison purposes.

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
