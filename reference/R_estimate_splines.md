# Estimate Rt via smoothing splines

This option estimates the effective reproduction number Rt using
penalized cubic basis splines.

## Usage

``` r
R_estimate_splines(
  knot_distance_global = 4 * 7,
  knot_distance_local = 7,
  R_start_prior_mu = 1,
  R_start_prior_sigma = 0.8,
  R_sd_local_prior_mu = 0,
  R_sd_local_prior_sd = 0.05,
  R_sd_global_prior_shape = 1,
  R_sd_global_prior_scale = 0.01,
  R_sd_global_change_distance = knot_distance_global,
  link = "inv_softplus",
  R_max = 6,
  modeldata = modeldata_init()
)
```

## Arguments

- knot_distance_global:

  Distance between spline breakpoints (knots) for the *global* spline in
  days. `EpiSewer` uses an ensemble of two splines to model Rt. The
  global spline models larger, long-term changes. Default is 4\*7, i.e.
  one knot every four week.

- knot_distance_local:

  Distance between spline breakpoints (knots) for the *local* spline in
  days. `EpiSewer` uses an ensemble of two splines to model Rt. The
  local spline models smaller, short-term changes. Default is 7, i.e.
  one knot each week.

- R_start_prior_mu:

  Prior (mean) on the initial reproduction number (intercept).

- R_start_prior_sigma:

  Prior (standard deviation) on the initial reproduction number
  (intercept).

- R_sd_local_prior_mu:

  Prior (mean) for the variation of the local (i.e. short-term) spline.
  This controls the standard deviation of a random walk over the
  coefficients of the local spline. The prior refers to the *daily*
  standard deviation and is thus independent of the knot distance.

- R_sd_local_prior_sd:

  Prior (standard deviation) for the variation of the local (i.e.
  short-term) spline. See `R_sd_local_prior_mu` for details.

- R_sd_global_prior_shape:

  Exponential-Gamma prior (shape) for the variation of the global (i.e.
  long-term) spline. This controls the standard deviation of a random
  walk over the coefficients of the global spline. The prior refers to
  the *daily* standard deviation and is thus independent of the knot
  distance. The exponential-Gamma prior is sparse, i.e. it has a strong
  peak towards zero and a long tail. Smaller shape parameters will lead
  to more sparseness, i.e. a longer tail. Note that when adjusting the
  shape, you will likely also have to adjust the scale. The variation of
  the global splines follows a change point model to allow for adaptive
  changes, see details for more explanation.

- R_sd_global_prior_scale:

  Exponential-Gamma prior (scale) for the variation of the global (i.e.
  long-term) spline. Larger scales will lead to more variability. See
  `R_sd_global_prior_shape` and the details for more explanation about
  this prior.

- R_sd_global_change_distance:

  Distance between changepoints used to model global variation in Rt.
  `EpiSewer` uses an adaptive model for the variation of the global
  spline, to model both time periods with stable and with volatile
  transmission dynamics. The default change point distance is equal to
  `knot_distance_global`. Making this smaller than the global knot
  distance will not have a large effect, however making this longer
  reduces the adaptability of the Rt variation. If set to zero, no
  change points are modeled, meaning zero adaptability.

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

`EpiSewer` uses a combination of two splines (global, for
larger/long-term changes and local, for smaller/short-term changes),
with a random walk on the coefficients of each spline. This allows a
highly adaptive yet regularized Rt model.

- The prior `R_start_prior` should reflect your expectation of Rt at the
  beginning of the time series. If your time series starts in the midst
  of an outbreak, you might want to use a prior with mean larger than 1
  for the intercept.

- The prior `R_sd_local_prior` on the standard deviation of the local
  spline coefficients should be interpreted in terms of daily additive
  changes (this is accurate around Rt=1, and becomes less accurate as Rt
  approaches 0 or its upper bound as defined by the `link` function).
  For example, a baseline half-normal prior with sd=0.05 allows a daily
  standard deviation between 0 and 0.1. A daily standard deviation of
  0.1 in turn roughly allows the spline coefficients to change by ±0.2
  (using the 2 sigma rule) each day. The daily standard deviation is
  summed up over the days between two knots to get the actual standard
  deviation of the coefficients. This way, the prior is independent of
  the chosen `knot_distance`. For example, if `knot_distance` is 7 days,
  and a constant daily standard deviation of 0.1 is estimated, the
  coefficients of two adjacent splines can differ by up to
  `0.2*sqrt(knot_distance)`, i.e. ±0.5. Note however that this does not
  directly translate into a change of Rt by ±0.5, as Rt is always the
  weighted sum of several basis functions at any given point. It will
  therefore change more gradually, depending on the distances between
  knots.

- The prior `R_sd_global_prior` on the standard deviation of the global
  spline coefficients is interpreted in the same way as the local
  spline, i.e. it describes daily additive changes. However, we us a
  *sparse* prior for the global variation, combined with a changepoint
  model. This yields a global spline that expects small changes most of
  the time but allows for occasional large changes. For example, during
  the height of an epidemic wave, countermeasures may lead to much
  faster changes in Rt than observable at other times. These differences
  in variability are accounted for using change points placed at regular
  intervals, with the daily global standard deviation evolving linearly
  between the change points. The values at the changepoints are modeled
  as independently distributed and follow the Exponential-Gamma (EG)
  prior distribution defined by `R_sd_global_prior`. The EG distribution
  is also known as Lomax distribution and corresponds to an exponential
  distribution with a Gamma distributed rate parameter.

The smoothness of the Rt estimates is influenced by a combination of the
global and local knot distances and by the priors on the global and
local standard deviation of the random walk on spline coefficients.
Placing knots further apart increases the smoothness of Rt estimates and
can speed up model fitting. The Rt time series remains surprisingly
flexible even at larger knot distances, but placing knots too far apart
can lead to inaccurate estimates. Note that the priors for the variation
of the global and local random walk also influence the uncertainty of Rt
estimates towards the present / date of estimation, when limited data
signal is available. Absent sufficient data signal, Rt estimates will
tend to stay at the current level (which corresponds to assuming
unchanged transmission dynamics).

The priors of this component have the following functional form:

- initial Rt (intercept): `Normal`

- daily standard deviation of the random walk over local spline
  coefficients: `Half-normal`

- daily standard deviation of the random walk over global spline
  coefficients: `Exponential-Gamma`

## See also

Other Rt models:
[`R_estimate_approx()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_approx.md),
[`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md),
[`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md),
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md),
[`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md),
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md),
[`R_estimate_smooth_derivative()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_smooth_derivative.md)
