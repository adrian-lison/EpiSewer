# Assume a limit of detection

Pathogen concentrations below a certain threshold may not be detectable
and thus erroneously measured as 0. This option adjusts for a known
limit of detection and includes zero measurements in the likelihood.

## Usage

``` r
LOD_assume(
  limit = NULL,
  prob = 0.95,
  LOD_type = "exponential",
  drop_prob = 1e-10,
  modeldata = modeldata_init()
)
```

## Arguments

- limit:

  Limit of detection. The concentration below which the pathogen cannot
  be detected with sufficient probability, i.e. the measurement may be
  zero although the pathogen is present in the sample.

- prob:

  What desired probability of detection does the limit refer to? Default
  is 95% (0.95): This means that the provided `limit` is the smallest
  concentration at which the pathogen can still be detected with over
  95% probability.

- LOD_type:

  The type of LOD model used. Currently, only "exponential" is
  supported. This models an exponentially decreasing probability of zero
  measurements / non-detection as a function of concentration. The
  exponential model can be derived from the statistical properties of
  dPCR, but should also work well for other quantification methods such
  as qPCR.

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

The limit of detection is modeled using a hurdle model. In effect, zero
measurements provide a signal that the concentration in the respective
sample was likely below the limit of detection, but we don't know what
the exact concentration was.

The limit of detection is specific to the quantification approach and
protocol. It is usually established from a dedicated lab experiment
(serial dilution experiment). It his here assumed that this experiment
did not cover a large fraction of the preprocessing noise to find an
optimal configuration for the exponential model.

If used together with
[`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md),
EpiSewer will also model the effect of pre-PCR noise on the LOD. This
means that the modeled LOD could be slightly higher than specified under
`limit`, depending on the estimated pre-PCR noise.

## See also

Visualize the assumed LOD as a function of concentration:
[`plot_LOD()`](https://adrian-lison.github.io/EpiSewer/reference/plot_LOD.md)

Other LOD models:
[`LOD_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md),
[`LOD_none()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_none.md)
