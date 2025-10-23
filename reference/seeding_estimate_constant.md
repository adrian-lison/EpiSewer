# Estimate constant seeding infections

This option estimates a constant number of infections at the start of
the modeled time period.

## Usage

``` r
seeding_estimate_constant(
  intercept_prior_q5 = NULL,
  intercept_prior_q95 = NULL,
  extend = TRUE,
  modeldata = modeldata_init()
)
```

## Arguments

- intercept_prior_q5:

  Prior (5% quantile) on the initial number of infections. Can be
  interpreted as an approximate lower bound. If NULL (default), this is
  computed from a crude empirical estimate of the number of cases (see
  details).

- intercept_prior_q95:

  Prior (95% quantile) on the initial number of infections. Can be
  interpreted as an approximate upper bound. If NULL (default), this is
  computed from a crude empirical estimate of the number of cases (see
  details).

- extend:

  Should the seeding phase be extended when concentrations are very low
  at the start of the measurement time series? If `TRUE`, then the
  seeding phase will be extended to the first date with three
  consecutive detects (i.e. non-zero measurements). Before that date,
  the reproduction number will be retrospectively computed based on the
  seeded infections. Explicit modeling of the Rt time series will only
  begin after the seeding phase. This option avoids sampling problems
  when estimating Rt from very low infection numbers during a period
  with many non-detects. Note that estimated reproduction numbers are
  not necessarily meaningful during periods with very low infection
  numbers, as transmission dynamics may be dominated by chance events
  and importations.

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

The seeding phase has the length of the maximum generation time (during
this time, the renewal model cannot be applied). It is here assumed that
the expected number of new infections stays constant over this time
period. This assumption can however be violated: While traditionally,
seeding refers to the first few (potentially imported) infections of an
epidemic, depending on what time period the model is fitted to, it may
also cover a different phase with stronger growth dynamics. Thus, if
your data starts in the middle of an epidemic wave, it is recommended to
use
[`seeding_estimate_growth()`](https://adrian-lison.github.io/EpiSewer/reference/seeding_estimate_growth.md)
instead of `seeding_estimate_constant()`.

If `intercept_prior_q5` or `intercept_prior_q95` are not specified by
the user, `EpiSewer` will compute a rough median empirical estimate of
the number of cases using the supplied wastewater measurements and
shedding assumptions, and then infer the missing quantiles based on
this. If none of the quantiles are provided, they are set to be roughly
1/10 and 10 times the empirical median estimate. We note that this is a
violation of Bayesian principles (data must not be used to inform
priors) - but a neglectable one, since it only ensures that the seeding
is modeled on the right order of magnitude and does not have relevant
impacts on later Rt estimates.

The priors of this component have the following functional form:

- initial number of infections (log scale): `Normal`
