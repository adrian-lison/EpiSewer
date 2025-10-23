# Estimate Rt via exponential smoothing

This option estimates the effective reproduction number over time using
exponential smoothing. It implements Holt's linear trend method with
damping through an innovations state space model with a level, trend,
and damping component.

## Usage

``` r
R_estimate_ets(
  R_start_prior_mu = 1,
  R_start_prior_sigma = 0.8,
  trend_prior_mu = 0,
  trend_prior_sigma = 0.1,
  sd_base_prior_mu = 0,
  sd_base_prior_sd = 0.025,
  sd_change_prior_shape = 0.5,
  sd_change_prior_scale = 1e-04,
  sd_change_distance = 7 * 26,
  link = "inv_softplus",
  R_max = 6,
  smooth_prior_mu = 0.5,
  smooth_prior_sigma = 0.05,
  trend_smooth_prior_mu = 0.5,
  trend_smooth_prior_sigma = 0.05,
  dampen_prior_mu = 0.9,
  dampen_prior_sigma = 0,
  differenced = FALSE,
  noncentered = TRUE,
  modeldata = modeldata_init()
)
```

## Arguments

- R_start_prior_mu:

  Prior (mean) on the initial value of Rt (level).

- R_start_prior_sigma:

  Prior (standard deviation) on the initial value of Rt (level).

- trend_prior_mu:

  Prior (mean) on the initial trend of Rt.

- trend_prior_sigma:

  Prior (standard deviation) on the initial trend of Rt.

- sd_base_prior_mu:

  Prior (mu) on the baseline standard deviation of the innovations.
  Please note that for consistency, the overall standard deviation of
  innovations will always be the baseline plus an additive component
  from `sd_change_prior` - even if no changepoints are modeled (see
  below).

- sd_base_prior_sd:

  Prior (standard deviation) on the baseline standard deviation of the
  innovations. See `sd_base_prior_mu` for details.

- sd_change_prior_shape:

  Exponential-Gamma prior (shape) on standard deviation additional to
  baseline. This prior describes the distribution of the standard
  deviation of Rt over time. EpiSewer will estimate a baseline standard
  deviation (see `sd_base_prior_mu`), and model additional variation on
  top of the baseline using a changepoint model. Please see the details
  for more explanation.

- sd_change_prior_scale:

  Exponential-Gamma prior (scale) on standard deviation additional to
  baseline. See `sd_change_prior_shape` and the details for more
  explanation.

- sd_change_distance:

  Distance between changepoints used to model additional variation in
  Rt. The default change point distance is 4 weeks. Very short
  changepoint distances must be chosen with care, as they can make the
  Rt time series too flexible. If set to zero, no change points are
  modeled.

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

- smooth_prior_mu:

  Prior (mean) on the smoothing parameter. Must be between 0 and 1.

- smooth_prior_sigma:

  Prior (standard deviation) on the smoothing parameter. If this is set
  to zero, the smoothing parameter will be fixed to `smooth_prior_mu`
  and not estimated. If positive, a beta prior with the corresponding
  mean and standard deviation is used.

- trend_smooth_prior_mu:

  Prior (mean) on the trend smoothing parameter. Must be between 0 and
  1.

- trend_smooth_prior_sigma:

  Prior (standard deviation) on the trend smoothing parameter. If this
  is set to zero, the trend smoothing parameter will be fixed to
  `trend_smooth_prior_mu` and not estimated. If positive, a beta prior
  with the corresponding mean and standard deviation is used.

- dampen_prior_mu:

  Prior (mean) on the damping parameter. Must be between 0 and 1.

- dampen_prior_sigma:

  Prior (standard deviation) on the damping parameter. If this is set to
  zero, the damping parameter will be fixed to `dampen_prior_mu` and not
  estimated. If positive, a beta prior with the corresponding mean and
  standard deviation is used.

- differenced:

  If `FALSE` (default), exponential smoothing is applied to the absolute
  Rt time series. If `TRUE`, it is instead applied to the differenced
  time series. This makes the level become the trend, and the trend
  become the curvature.

- noncentered:

  If `TRUE` (default), a non-centered parameterization is used to model
  the innovations in the state space process (for better sampling
  efficiency).

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

The innovations state space model consists of three components: a level,
a trend, and a damping component.

- The level is smoothed based on the levels from earlier time steps,
  with exponentially decaying weights, as controlled by a smoothing
  parameter (often called alpha). Note that *smaller* values of `alpha`
  indicate *stronger* smoothing. In particular, `alpha = 1` means that
  only the last level is used.

- The trend is smoothed based on the trends from earlier time steps,
  with exponentially decaying weights, as controlled by a trend
  smoothing parameter (often called beta). Note that *smaller* values of
  `beta` indicate *stronger* smoothing. In particular, `beta = 1` means
  that only the last trend is used.

- The damping determines how long a previous trend continues into the
  future before it levels of to a stationary time series. The strength
  of damping is controlled by a damping parameter (often called phi).
  Note that *smaller* values of `phi` indicate *stronger* damping. In
  particular, `phi = 1` means no damping. Values below `phi = 0.8` are
  seldom in practice as the damping becomes very strong.

Often, `alpha`, `beta`, and `phi` are jointly unidentifiable. It may
therefore be necessary to fix at least one of the parameters (typically
`phi`) or supply strong priors.

Note that the smoothness of retrospective Rt estimates is often more
influenced by the prior on the standard deviation of innovations than
the smoothing and trend smoothing parameters. The smoothing parameters
mostly have an influence on estimates close to the present / date of
estimation, when limited data signal is available. Here, the standard
deviation of the innovations influences how uncertain Rt estimates are
close to the present.

The variability of Rt can change over time. For example, during the
height of an epidemic wave, countermeasures may lead to much faster
changes in Rt than observable at other times (baseline). This potential
additional variability is accounted for using change points placed at
regular intervals. The additional standard deviation of the state space
model innovations on top of the baseline then evolves linearly between
the change points. The additional variation defined at the changepoints
is modeled as independently distributed and following a Lomax
distribution, also known as Exponential-Gamma (EG) distribution. This is
an exponential distribution where the rate is Gamma distributed. The
prior `sd_change_prior` defines the shape and scale of this Gamma
distribution. The distribution has a strong peak towards zero and a long
tail. This regularizes the estimated deviations from the baseline
standard deviation - most deviations are small, but during special time
periods, the deviation might also be larger.

The priors of this component have the following functional form:

- initial level of Rt: `Normal`

- initial trend of Rt: `Normal`

- baseline standard deviation of innovations: `Half-normal`

- additional standard deviation at changepoints: `Exponential-Gamma`

- smoothing parameter: `Beta`

- trend smoothing parameter: `Beta`

- damping parameter: `Beta`

## See also

Other Rt models:
[`R_estimate_approx()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_approx.md),
[`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md),
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md),
[`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md),
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md),
[`R_estimate_smooth_derivative()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_smooth_derivative.md),
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
