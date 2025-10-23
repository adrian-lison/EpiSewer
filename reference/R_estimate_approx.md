# Estimate Rt using an approximation of the generative renewal model

**\[deprecated\]**

This function is deprecated. Please use
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md)
or other fully generative methods instead.

This option estimates the effective reproduction number Rt using an
approximation of the generative renewal model. It can be considerably
faster than the fully generative options, but is often less exact. It is
recommended to check the results of this method against a fully
generative version like
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md)
before using it for real-time R estimation. See the details for an
explanation of the method and its limitations.

## Usage

``` r
R_estimate_approx(
  inf_sd_prior_mu = 0.05,
  inf_sd_prior_sigma = 0.025,
  inf_smooth = 0.75,
  inf_trend_smooth = 0.75,
  inf_trend_dampen = 0.9,
  knot_distance = 7,
  spline_degree = 3,
  R_window = 1,
  modeldata = modeldata_init()
)
```

## Arguments

- inf_sd_prior_mu:

  Prior (mean) on the standard deviation of daily changes in infections.
  Infection numbers are modeled on the log scale, thus the standard
  deviation can be roughly interpreted in terms of percentage changes. A
  higher standard deviation will lead to more volatile infections and
  more uncertainty in R estimates. The default
  (`inf_sd_prior_mu = 0.05`) allows daily growth/decline rates of up
  +-10%. We suggest to only decrease this if you are very certain that
  the modeled scenario has lower maximum growth/decline rates, and to
  increase this if you expect a more extreme scenario. Note that the
  smoothness of the infection curve is also influenced by the
  exponential smoothing parameters (`inf_smooth`, `inf_trend_smooth`).
  Still, `inf_sd_prior_mu` is the primary hyperparameter to adjust if
  you think the estimated infection trajectory is too volatile / not
  flexible enough.

- inf_sd_prior_sigma:

  Prior (standard deviation) on standard deviation of daily changes of
  infections. This reflects uncertainty around `inf_sd_prior_mu`. The
  default (`inf_sd_prior_mu=0.025`) allows the standard deviation to be
  roughly 0.5 higher or lower than `inf_sd_prior_mu`. For the default
  setting with `inf_sd_prior_mu = 0.05`, this means that growth/decline
  rates of up to +-20% are still supported by the prior.

- inf_smooth:

  Smoothing parameter (alpha) for infections. The infection time series
  is smoothed based on earlier time steps, with exponentially decaying
  weights, as controlled by the alpha parameter. The default is 0.5,
  *smaller* values will lead to *stronger* smoothing, and *larger*
  values to *less* smoothing.

- inf_trend_smooth:

  Trend smoothing parameter (beta) for infections. The exponential trend
  in infections is also smoothed based on the trend over earlier time
  steps, with exponentially decaying weights, as controlled by the beta
  parameter. Default is 0.5, *smaller* values will lead to *stronger*
  smoothing of the trend (useful for longer generation times), and
  *larger* values to *less* smoothing of the trend (useful for shorter
  generation times).

- inf_trend_dampen:

  Trend damping parameter (phi) for infections. The exponential trend in
  infections is dampened over time (exponential growth will not continue
  for ever). This mainly influences the trend towards the present. For
  the default value (0.9), the trend will only roughly half as strong
  after two weeks. Note that increasing phi to close to 1 can reduce
  sampling speed.

- knot_distance:

  The infection time series is smoothed using penalized B-splines. The
  distance between spline breakpoints (knots) is 7 days by default.
  Placing knots further apart increases the smoothness of the infection
  time series and can speed up model fitting. Placing knots closer to
  each other increases the volatility of the infections. Note however
  that the uncertainty of R estimates produced by `R_estimate_approx()`
  is also sensitive to the smoothness of the infection time series, i.e.
  changing the knot distance can make the R estimates appear more or
  less uncertain. This effect is however rather small since the
  regularization of the spline coefficients adjusts for the knot
  distance.

- spline_degree:

  Degree of the spline polynomials (default is 3 for cubic splines).

- R_window:

  Smoothing window size for R estimates. Default is 1 (i.e. no
  smoothing). This is provided for compatibility with the EpiEstim
  method by Cori et al., which assumes that R stays constant over a
  certain time window. However, this can artificially reduce the
  uncertainty of the R estimates and is therefore generally not
  recommended. It also means that R cannot be estimated up to the
  present because the smoothing window must be centered.

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

This method uses an approximation of the renewal model to estimate the
effective reproduction number. Instead of generating the infections
through a renewal process, the infection time series is estimated using
a non-parametric smoothing prior that mimics the characteristics of a
renewal process. Rt is then computed from the infection time series by
applying the classical renewal equation. This still gives Bayesian
estimates of Rt and the infection time series, but is less exact than a
fully generative model.

To smooth the infection time series, penalized B-splines similar to the
ones in
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
are used, but they are applied to the expected number of infections on
the log scale, not Rt. The spline coefficients are themselves smoothed
using an exponential smoothing prior with a trend component. The trend
component is important as it reflects the default assumption of constant
transmission dynamics that leads to exponential growth in infections.
This ensures consistent results of `R_estimate_approx()` also towards
the present.

The method can also model infection noise as specified by
[`infection_noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/infection_noise_estimate.md).
To account for the autocorrelation of infection noise, a correction that
applies the renewal model specifically to the noise component is used.
This correction is most accurate when R is close to 1 and less accurate
when it is far from 1.

The main limitation of `R_estimate_approx()` is that its priors and
hyperparameters are less interpretable and thus difficult to specify. In
particular the assumed generation time distribution is only used for R
estimation but does not inform the smoothness of the expected infection
time series. Instead, the smoothness is determined by the spline
settings and exponential smoothing parameters, which are difficult to
translate into a generation time assumption. The best approach is
therefore to compare estimates from `R_estimate_approx()` with those
from another method like
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
for a pathogen of interest and to carefully adjust the smoothing
hyperparameters if needed.

## See also

Other Rt models:
[`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md),
[`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md),
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md),
[`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md),
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md),
[`R_estimate_smooth_derivative()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_smooth_derivative.md),
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
