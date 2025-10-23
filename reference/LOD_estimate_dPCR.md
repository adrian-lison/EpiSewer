# Estimate a limit of detection model for digital PCR data

Pathogen concentrations below a certain threshold may not be detectable
and thus erroneously measured as 0. This option adjusts for a limit of
detection based on the statistical properties of digital PCR (dPCR) and
includes zero measurements in the likelihood.

In effect, zero measurements provide a signal that the concentration in
the respective sample was likely below the limit of detection, but we
don't know what the exact concentration was.

## Usage

``` r
LOD_estimate_dPCR(drop_prob = 1e-10, modeldata = modeldata_init())
```

## Arguments

- drop_prob:

  Probability for non-detection below which likelihood contributions of
  observed concentrations are dropped from LOD model. This avoids
  numerical issues of the LOD model at high concentrations (very small
  non-detection probabilities) that can otherwise affect sampling speed.
  Since these likelihood contributions will be virtually zero for almost
  all samples anyway, parameter estimates are practically not affected.

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

The limit of detection is modeled using a hurdle model. The model uses
the number of partitions in the dPCR reaction and the conversion factor
as defined and estimated by
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md).
It can therefore only be used together with
`noise = noise_estimate_dPCR()` in
[`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md).

## See also

Other LOD models:
[`LOD_assume()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_assume.md),
[`LOD_none()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_none.md)
