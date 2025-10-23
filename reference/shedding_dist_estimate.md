# Estimate an uncertain shedding load distribution

This option models a parametric shedding load distribution. It is useful
for representing uncertainty about the mean and variation of the
shedding profile.

## Usage

``` r
shedding_dist_estimate(
  shedding_dist_mean_prior_mean = NULL,
  shedding_dist_mean_prior_sd = NULL,
  shedding_dist_cv_prior_mean = NULL,
  shedding_dist_cv_prior_sd = NULL,
  shedding_dist_type = "gamma",
  shedding_reference = NULL,
  prior_weights = NULL,
  weight_alpha = 1,
  modeldata = modeldata_init()
)
```

## Arguments

- shedding_dist_mean_prior_mean:

  Prior (mean) for the mean of the shedding load distribution. Can also
  be a vector to represent several priors, that will be mixed using a
  Dirichlet prior. This is useful when combining several shedding load
  distributions from different studies or sensitivity analyses.

- shedding_dist_mean_prior_sd:

  Prior (standard deviation) for the mean of the shedding load
  distribution. Can also be a vector, see
  `shedding_dist_mean_prior_mean`.

- shedding_dist_cv_prior_mean:

  Prior (mean) for the coefficient of variation (i.e. standard deviation
  relative to the mean) of the shedding load distribution. Can also be a
  vector, see `shedding_dist_mean_prior_mean`.

- shedding_dist_cv_prior_sd:

  Prior (standard deviation) for the coefficient of variation of the
  shedding load distribution. Can also be a vector, see
  `shedding_dist_mean_prior_mean`.

- shedding_dist_type:

  The parametric distribution that should be modeled. Supported are
  "gamma", "exponential", and "lognormal".

- shedding_reference:

  Is the shedding load distribution relative to the day of `"infection"`
  or the day of `"symptom_onset"`? This is important because shedding
  load distributions provided in the literature are sometimes by days
  since infection and sometimes by days since symptom onset. If
  `shedding_reference="symptom_onset"`, EpiSewer also needs information
  about the incubation period distribution (see
  [`incubation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/incubation_dist_assume.md)).

- prior_weights:

  A numeric vector of the same length as the priors for the mean and cv,
  with weights for the different priors. If `NULL`, the priors are given
  equal weight. Will be normalized to sum to the number of
  distributions.

- weight_alpha:

  Concentration parameter of the Dirichlet prior for the mixture
  probabilities. The default is `weight_alpha=1`, i.e. a uniform
  distribution over all combinations of weights. This will sample
  various mixtures of the mean and cv prior distributions provided. If a
  smaller `weight_alpha` is chosen, the prior distributions are rather
  considered separately, i.e. in each posterior sample, one of the prior
  distributions is given almost all the weight.

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
