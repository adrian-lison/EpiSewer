# Visually compare prior and posterior of a model parameter

Visually compare prior and posterior of a model parameter

## Usage

``` r
plot_prior_posterior(result, param_name)
```

## Arguments

- result:

  Results object returned by
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  after model fitting. In contrast to other plotting functions, this
  cannot be a list of multiple result objects, because prior-posterior
  plots of multiple results simultaneously are currently not supported.

- param_name:

  Name of the single parameter to be plotted. This can either be the raw
  name of the parameter in the stan model, or it can be the semantic
  name from the single parameter dictionary (see details). Note that
  this must be a single parameter, it cannot be a parameter array. Also,
  only original parameters (not transformed parameters) for which a
  prior can be specified in `EpiSewer` are supported.

## Value

A plot showing the density of the prior (grey) and posterior (blue) for
the respective parameter. Can be further manipulated using `ggplot2`
functions to adjust themes and scales, and to add further geoms.

## Details

The following parameters can be visualized (if in the model):

- `measurement_noise_cv` (nu_upsilon_a): Coefficient of variation
  (measurement noise)

- `dPCR_total_partitions` (nu_upsilon_b_mu): Average total number of
  partitions in dPCR

- `dPCR_partition_variation` (nu_upsilon_b_cv): Partition number
  variation in dPCR

- `dPCR_conversion_factor` (nu_upsilon_c): Conversion factor in dPCR

- `pre_replicate_cv` (nu_psi): Coefficient of variation (pre-PCR noise)

- `load_variation_cv` (nu_zeta): Individual-level coefficient of load
  variation

- `infection_overdispersion` (I_xi): Overdispersion of infections

- `seeding_intercept` (iota_log_seed_intercept): Initial number of
  infections
