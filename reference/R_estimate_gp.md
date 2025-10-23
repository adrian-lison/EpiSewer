# Estimate Rt using Gaussian processes

This option estimates the effective reproduction number Rt over time
using a Gaussian process (GP) model. There are two GPs: one for the
long-term trend in Rt, and one for short-term deviations from this
trend.

## Usage

``` r
R_estimate_gp(
  R_intercept_prior_mu = 1,
  R_intercept_prior_sigma = 0,
  length_scale_prior_mu = 7 * 3,
  length_scale_prior_sigma = 7/2,
  magnitude_prior_mu = 0.2,
  magnitude_prior_sigma = 0.05,
  long_length_scale_prior_mu = 7 * 4 * 3,
  long_length_scale_prior_sigma = 7,
  long_magnitude_prior_mu = 0.4,
  long_magnitude_prior_sigma = 0.1,
  matern_nu = c(3/2, 5/2, 1/2),
  boundary_factor = 3,
  n_basis_factor = 3.42,
  link = "inv_softplus",
  R_max = 6,
  modeldata = modeldata_init()
)
```

## Arguments

- R_intercept_prior_mu:

  Prior (mean) for the intercept of Rt. Should be set to 1 unless you
  have a clear a priori expectation of the average Rt during the modeled
  time period.

- R_intercept_prior_sigma:

  Prior (standard deviation) for the intercept of Rt. By default, we fix
  `R_intercept` to 1 by setting `R_intercept_prior_sigma=0`. This is
  because deviations from Rt=1 are already captured by the long-term
  trend component (see below).

- length_scale_prior_mu:

  Prior (mean) on the length scale of the short-term Gaussian process
  (in days). This influences the smoothness of Rt. A higher length scale
  means that the Rt will change more slowly. Choosing a length scale
  that is too short can lead to overfitting of the Rt trajectory, while
  choosing a length scale that is too long can result in unrealistically
  smooth Rt estimates.

- length_scale_prior_sigma:

  Prior (standard deviation) on the length scale of the short-term
  Gaussian process. Set to zero to fix the length scale.

- magnitude_prior_mu:

  Prior (mean) on the magnitude of the short-term Gaussian process. Can
  be approximately interpreted as the marginal standard deviation of Rt.
  A higher magnitude allows more extreme Rt values.

- magnitude_prior_sigma:

  Prior (standard deviation) on the magnitude of the short-term Gaussian
  process. Set to zero to fix the magnitude.

- long_length_scale_prior_mu:

  Prior (mean) on the length scale of the long-term Gaussian process (in
  days). This should be quite long, at least several times the mean
  shedding delay of the pathogen, and significantly larger than the mean
  prior for the short-term GP (`length_scale_prior_mu`).

- long_length_scale_prior_sigma:

  Prior (standard deviation) on the length scale of the long-term
  Gaussian process. Set to zero to fix the length scale.

- long_magnitude_prior_mu:

  Prior (mean) on the magnitude of the long-term Gaussian process. Can
  be approximately interpreted as the marginal standard deviation of the
  long-term Rt trend. A higher magnitude allows more extreme Rt values.

- long_magnitude_prior_sigma:

  Prior (standard deviation) on the magnitude of the long-term Gaussian
  process. Set to zero to fix the magnitude.

- matern_nu:

  The smoothness parameter of the Matern kernel. The default is 3/2,
  other possible choices are 5/2 (more smooth) and 1/2 (less smooth).
  However, we recommend tuning smoothness primarily using the length
  scale priors.

- boundary_factor:

  The boundary factor used in the Gaussian process approximation. The
  default (`boundary_factor = 3`) is higher than the minimum
  recommendation from Riutort-Mayol et al. to ensure accurate real-time
  estimation. Lower values can lead to boundary effects close to the
  present and start of the time series and are thus not recommended.
  When modeling very long shedding delays, the boundary factor might
  need to be further increased. Note that this will also automatically
  increase the number of basis functions used and thereby slow down
  sampling.

- n_basis_factor:

  Factor used to automatically determine `m`, the number of basis
  functions in the Gaussian process approximation. Based on
  recommendations from Riutort-Mayol et al., we let
  `m = n_basis_factor * c / (l / S)`, where `c` is the
  `boundary_factor`, `l` is the length scale of the GP (we use the 5%
  quantile of `gp_length_prior` for a conservative result), and `S` is
  the maximum absolute value of the zero-centered input space, which is
  `(n-1)/2` where `n` is the number of Rt time steps modeled. This means
  that the number of basis function automatically adapts to the length
  scale and number of observations in the model. For the default
  settings (`boundary_factor = 3` and `n_basis_factor = 3.42`) and a
  lower length scale of 3 weeks, this corresponds to approx. 0.25x the
  number of Rt time points modeled. Increasing `n_basis_factor` will
  make the approximation more accurate by proportionally increasing the
  number of basis functions, but can slow down sampling.

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

The estimated Rt trajectory is primarily influenced by the priors for
the *length scale* and the *magnitude* of the short-term and long-term
Gaussian processes. We recommend adjusting the *magnitude* when maximum
Rt values seem to be too low or too high. In contrast, we recommend
adjusting the *length scale* when the Rt trajectory seems to have too
little resolution (try shorter length scales) or is overfitting on noise
(try longer length scales). Note that, to ensure identifiability, the
prior mean for the length scale of the long-term GP should always be
significantly larger than the prior mean for the short-term GP.

The Gaussian process is modeled using a Hilbert space approximation as
described in Riutort-Mayol et al. (2023). This allows for fast inference
without significant loss of accuracy. See the reference for more detail
on choosing an adequate boundary factor and basis function factor. The
`EpiSewer` implementation of the approximate Gaussian process is
strongly inspired by the implementation in the
[EpiNow2](https://github.com/epiforecasts/EpiNow2/) package.

The priors of this component have the following functional form:

- R_intercept (intercept): `Normal`

- magnitude (magnitude): `Truncated Normal`

- length_scale (length scale): `Truncated Normal`

- long_magnitude (magnitude for long-term trend): `Truncated Normal`

- long_length_scale (length scale for long-term trend):
  `Truncated Normal`

## References

Riutort-Mayol, G., Bürkner, PC., Andersen, M.R. et al. Practical Hilbert
space approximate Bayesian Gaussian processes for probabilistic
programming. *Stat Comput* 33, 17 (2023).
<https://doi.org/10.1007/s11222-022-10167-2>

Abbott S, Hellewell J, Thompson RN et al. Estimating the time-varying
reproduction number of SARS-CoV-2 using national and subnational case
counts. *Wellcome Open Res* 5, 112 (2020).
<https://doi.org/10.12688/wellcomeopenres.16006.2>

## See also

Other Rt models:
[`R_estimate_approx()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_approx.md),
[`R_estimate_changepoint_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_changepoint_splines.md),
[`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md),
[`R_estimate_piecewise()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_piecewise.md),
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md),
[`R_estimate_smooth_derivative()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_smooth_derivative.md),
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
