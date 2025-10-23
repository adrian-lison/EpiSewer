# Estimate Rt via a changepoint spline model

This option models the effective reproduction number Rt over time using
cubic splines that are regularized via a linear changepoint model. The
Rt trajectory will thus follow a smoothed linear trend, with time points
of trend changes estimated from the data using an approximate
changepoint model. This approach offers high flexibility of the Rt
trajectory while avoiding overfitting on noise.

## Usage

``` r
R_estimate_changepoint_splines(
  R_start_prior_mu = 1,
  R_start_prior_sigma = 0.8,
  changepoint_max_distance = 3 * 5,
  changepoint_min_distance = 3 * 2,
  trend_prior_shape = 50,
  trend_prior_scale = 1,
  trend_change_tolerance = 0.01,
  spline_knot_distance = 3,
  link = "inv_softplus",
  R_max = 6,
  strictness_alpha = 0.5,
  modeldata = modeldata_init()
)
```

## Arguments

- R_start_prior_mu:

  Prior (mean) on the initial reproduction number (intercept).

- R_start_prior_sigma:

  Prior (standard deviation) on the initial reproduction number
  (intercept).

- changepoint_max_distance:

  Maximum distance (in days) between changes of the Rt trend. This
  setting guarantees that two consecutive changes in the Rt trend can be
  captured if they are `changepoint_max_distance` days or more apart.
  Faster changes are not guaranteed to be captured.

- changepoint_min_distance:

  Minimum distance (in days) between changes of the Rt trend. This
  setting guarantees that two consecutive changes in the Rt trend are at
  least `changepoint_min_distance` days apart. This avoids Rt
  trajectories that are unrealistically volatile. For example, a minimum
  distance of 7 implies that you don't expect Rt to significantly change
  its trend more than once within a week.

- trend_prior_shape:

  Exponential-Gamma (EG) prior (shape) for the trend in Rt. At each
  estimated changepoint, the Rt trend is sampled from a normal
  distribution with standard deviation given by the EG prior. This prior
  has a strong peak towards zero and a long tail. In other words, while
  we expect Rt to remain stable most of the time, this prior also allows
  for occasional strong trends. Smaller shape parameters will lead to a
  longer tail, hence more extreme trends are supported. Note that when
  adjusting the shape, you will likely also have to adjust the scale.
  See details for more advice on choosing a suitable prior.

- trend_prior_scale:

  Exponential-Gamma (EG) prior (scale) for the trend in Rt. See
  `change_prior_shape` above for an explanation. Larger scales will lead
  to more variability: a doubling of the scale roughly corresponds to a
  doubling of all quantiles of the prior.

- trend_change_tolerance:

  Tolerance for "negligible" trend changes. Differences in the trend
  that are smaller than `change_tolerance` are ignored by
  `changepoint_min_distance`, i.e. they can also occur closer to each
  other. This tolerance gives the model more flexibility in placing
  changepoints with large trend changes.

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

- strictness_alpha:

  The concentration parameter of the Dirichlet prior for the changepoint
  positions. Choosing smaller values of `strictness_alpha` will lead to
  more strict changepoints. Note that choosing small values of
  `strictness_alpha` can impede MCMC sampling.

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

The Exponential-Gamma (EG) prior on the Rt trend is parameterized via
the arguments `change_prior_shape` and `change_prior_scale`. It has a
long tail to support large trends while keeping the Rt variation low
most of the time. The default configuration should work well in most
contexts except for really extreme changes in Rt over a short time
window. To check the quantiles of your prior, you can use the function
[`qexpgamma()`](https://adrian-lison.github.io/EpiSewer/reference/qexpgamma.md)
with corresponding shape and scale parameters.

If you need to adjust the overall variation, you can adjust the
`change_prior_scale` parameter. A doubling of the scales roughly
corresponds to a doubling of the quantiles. For example, when the 95%
quantile is 0.2 for a given scale and you double that scale, the 95%
quantile will be at 0.4.

If you need to support more extreme changes, you can decrease the
`change_prior_shape` parameter, which will emphasize the long-tail
behavior of the prior. Note however that this will substantially
increase all quantiles of the prior, so you will also have to decrease
the scale parameter to achieve a similar level of day-to-day variation.

## See also

Other Rt models:
[`R_estimate_approx()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_approx.md),
[`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md),
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md),
[`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md),
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md),
[`R_estimate_smooth_derivative()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_smooth_derivative.md),
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
