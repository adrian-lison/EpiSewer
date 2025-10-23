# Estimate Rt with a smooth derivative

This option estimates the effective reproduction number Rt over time as
a smooth trend. It uses sparse smoothing splines placed on the
first-order differences of Rt.

## Usage

``` r
R_estimate_smooth_derivative(
  R_start_prior_mu = 1,
  R_start_prior_sigma = 0.8,
  spline_knot_distance = 14,
  trend_prior_shape = 5,
  trend_prior_scale = 0.005,
  link = "inv_softplus",
  R_max = 6,
  modeldata = modeldata_init()
)
```

## Arguments

- R_start_prior_mu:

  Prior (mean) on the initial reproduction number (intercept).

- R_start_prior_sigma:

  Prior (standard deviation) on the initial reproduction number
  (intercept).

- spline_knot_distance:

  Distance (in days) between spline knots for the penalized smoothing
  splines. Shorter distances increase flexibility of the Rt trajectory
  but also the daily Rt uncertainty.

- trend_prior_shape:

  Normal-Exponential-Gamma (NEG) prior (shape) for the Rt trend. The NEG
  prior is sparse, i.e. it has a strong peak at zero (transmission
  remains unchanged) and long tails (allows for occasional large
  positive or negative changes). Smaller shape parameters will lead to
  more sparseness (i.e. fewer but sharper changes in Rt). Note that when
  adjusting the shape, you will likely also have to adjust the scale.

- trend_prior_scale:

  Normal-Exponential-Gamma (NEG) prior (scale) for the strength of the
  Rt trend. See `sharp_changes_prior_shape` above for an explanation.
  Larger scales will lead to more Rt variability.

- link:

  Link function. Currently supported are `inv_softplus` (default) and
  `scaled_logit`. Both of these links are configured to behave
  approximately like the identity function around R=1, but become
  increasingly non-linear below (and in the case of `scaled_logit` also
  above) R=1.

- R_max:

  If `link=scaled_logit` is used, a maximum reproduction number must be
  assumed. This should be higher than any realistic R value for the
  modeled pathogen. Default is 6.

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

The priors of this component have the following functional form:

- R_start (intercept): `Normal`

- trend_prior: `Normal-Exponential-Gamma`

## See also

Other Rt models:
[`R_estimate_approx()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_approx.md),
[`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md),
[`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md),
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md),
[`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md),
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md),
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
