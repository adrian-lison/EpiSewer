Detailed EpiSewer example
================

This vignette repeats the introductory example from the
[README](../README.md), but shows how to specify the different modules
and components in `EpiSewer` explicitly.

💡 For a conceptual overview over modules and components, we recommend
to first read the [model specification](model-specification.md)
vignette.

### Loading the package

We assume that `EpiSewer` is already successfully installed (see
[README](../README.md#Installing%20the%20package)).

``` r
library(EpiSewer)
library(data.table)
library(ggplot2)
```

### Data

This is again the wastewater data are from Zurich, Switzerland, provided
by EAWAG (Swiss Federal Institute of Aquatic Science and Technology).

``` r
data_zurich <- SARS_CoV_2_Zurich
```

As in the introductory example, we keep only measurements that were made
on Mondays and Thursdays to test `EpiSewer` on sparse data.

``` r
measurements_sparse <- data_zurich$measurements[,weekday := weekdays(data_zurich$measurements$date)][weekday %in% c("Monday","Thursday"),]
head(measurements_sparse, 20)
#>           date concentration  weekday
#>  1: 2022-01-03      455.7580   Monday
#>  2: 2022-01-06      330.7298 Thursday
#>  3: 2022-01-10      387.6885   Monday
#>  4: 2022-01-13      791.1111 Thursday
#>  5: 2022-01-17      551.7701   Monday
#>  6: 2022-01-20      643.9910 Thursday
#>  7: 2022-01-24      741.9150   Monday
#>  8: 2022-01-27      770.1810 Thursday
#>  9: 2022-01-31      627.1725   Monday
#> 10: 2022-02-03      561.2913 Thursday
#> 11: 2022-02-07      357.1349   Monday
#> 12: 2022-02-10      540.7527 Thursday
#> 13: 2022-02-14            NA   Monday
#> 14: 2022-02-17      554.2492 Thursday
#> 15: 2022-02-21      414.7324   Monday
#> 16: 2022-02-24      784.3849 Thursday
#> 17: 2022-02-28      732.9672   Monday
#> 18: 2022-03-03     1376.6457 Thursday
#> 19: 2022-03-07     1420.4823   Monday
#> 20: 2022-03-10     2128.1925 Thursday
```

### Modeling

In the [README](../README.md), we used a shortcut to collect all
observation data via `sewer_data()` and all assumptions via
`sewer_assumptions()`. These are convenience functions to make the
`EpiSewer` command more concise. In this example however, we want to be
explicit and will supply the data and assumptions to the relevant
modules and components (see the [model
specification](model-specification.md) vignette).

#### Measurements

We start with the first module, the observed measurements. The module
can be specified as follows:

``` r
ww_measurements <- model_measurements(
  concentrations = concentrations_observe(measurements = measurements_sparse),
  noise = noise_estimate(),
  LOD = LOD_none()
)
```

As we see, the module is defined using the `model_measurements()`
function. The arguments of this function represent the different module
components. In this case, there is:

- a `concentrations` component: defines pathogen concentrations in the
  wastewater
- a `noise` component: defines measurement noise
- a `LOD` component: defines a limit of detection in wastewater samples

We used a specific modeling function for each of the components. The
suffix of each function explains what we modeled:

- `concentrations_observe`: the concentrations are modeled as
  observation data to which the model can be fitted
- `noise_estimate`: the degree of measurement noise is estimated from
  the data
- `LOD_none`: we do not model a limit of detection

We can also print the module object to see what we have modeled:

``` r
ww_measurements
#> measurements
#>  |- concentrations_observe
#>  |- noise_estimate
#>  |- LOD_none
```

❗ Note that the modeling functions above have further arguments, and
`EpiSewer` will use default settings if they are not specified. For
example, we could set a custom prior on the measurement noise as
follows:

``` r
ww_measurements2 <- model_measurements(
  concentrations = concentrations_observe(measurements = measurements_sparse),
  noise = noise_estimate(cv_prior_mu = 0.5, cv_prior_sigma = 0.1), # prior on coefficient of variation
  LOD = LOD_none()
)
```

➡️ To find out which options are available for a given modeling
function, you can consult its function documentation.

Before we move to the next module, remember that there can be multiple
modeling options for a component, and `EpiSewer` offers one modeling
function for each option. To see all available options, you can use
`component_functions()`:

``` r
component_functions("LOD")
#> [1] "LOD_none()"           "LOD_assume()"         "LOD_estimate_dPCR()"
```

For example, we could assume a limit of detection of `2.56 gc/mL` (based
on a lab experiment by EAWAG, this was the lowest concentration at which
still 95% of replicates were positive) as follows:

``` r
ww_measurements3 <- model_measurements(
  concentrations = concentrations_observe(measurements = measurements_sparse),
  noise = noise_estimate(cv_prior_mu = 0.5, cv_prior_sigma = 0.1), # prior on coefficient of variation
  LOD = LOD_assume(limit = 2.56, prob = 0.95)
)
```

#### Sampling

Next is the `sampling` module. It models the process of sampling from
the wastewater:

``` r
ww_sampling <- model_sampling(
  sample_effects = sample_effects_none()
)
```

The sampling module currently has only one component, which allows to
model the effects of sample characteristics (like time from sampling
until processing, which can lead to degradation) on the detectable
pathogen concentrations. In this simple example, we do not model such
effects.

#### Sewage

The `sewage` module describes aspects of the sewage system:

``` r
ww_sewage <- model_sewage(
  flows = flows_observe(flows = data_zurich$flows),
  residence_dist = residence_dist_assume(residence_dist = c(1))
)
```

It includes two components:

- `flows`: Daily flow volumes at the sampling site. Here we use the
  measured flows at the treatment plant in Zurich, specified using
  `flows_observe()`.
- `residence_dist`: Sewer residence time distribution for pathogen
  particles. In large sewage systems, particles may travel longer than a
  day depending on where and when they were shed into the wastewater.
  This component allows to model the corresponding delay distribution.
  We here use the modeling function `residence_dist_assume()` to
  indicate that we make an assumption about the residence time. In our
  case, we assume that particles arrive at the sampling site within the
  same day (indicated by the vector `c(1)`, which puts 100% probability
  on a residence time of 0 days).

#### Shedding

The `shedding` module describes how infected individuals shed pathogen
particles into the wastewater:

``` r
ww_shedding <- model_shedding(
  incubation_dist = incubation_dist_assume(get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4, maxX = 10)),
  shedding_dist = shedding_dist_assume(get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397, maxX = 30)),
  load_per_case = load_per_case_assume(6e+11),
  load_variation = load_variation_none()
)
```

The shedding module requires a number of assumptions (see explanations
in the [README](../README.md)):

- `incubation_dist` and `shedding_dist`: We use the same incubation
  period and shedding load distribution as before.
- `load_per_case`: We again assume an average shedding load of
  `6e+11 gc/person` as before.
- `load_variation`: `EpiSewer` offers the option to model
  individual-level variation in shedding. However, we keep it simple
  here and do not model such variation (`load_variation_none()`).

#### Infections

Last but not least, the `infections` module describes the underlying
process leading to infected individuals:

``` r
ww_infections <- model_infections(
  generation_dist = generation_dist_assume(get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12)),
  R = R_estimate_splines(),
  seeding = seeding_estimate_rw(),
  infection_noise = infection_noise_estimate()
)
```

The module has the following epidemiological components:

- `generation_dist`: We use the same generation time distribution as
  before.
- `R`: This the effective reproduction number, our main parameter of
  interest. We here estimate it via smoothing B-splines
  (`R_estimate_splines()`), but there are other smoothing approaches
  available (see `component_functions("R")`).
- `seeding`: The renewal model used by `EpiSewer` requires a seeding
  phase during which the reproduction number cannot be modeled. We thus
  let the initial infections be estimated using a simple random walk
  model.
- `infection_noise`: We also model noise in the infection process,
  meaning that `EpiSewer` will fit a stochastic infection model. This
  makes sense in most cases.

#### Side notes

We have now specified all five `EpiSewer` modules. While this seems like
a lot of code, note that the majority of component specifications are
optional and have sensible defaults. In practice, you only need to
explicitly write out components that need certain assumptions or that
you want to customize. In the present example, we also (unnecessarily)
specified the `EpiSewer` defaults. For example, the `infections` module
described above could simply be written as

``` r
ww_infections <- model_infections(
  generation_dist = generation_dist_assume(get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12))
)
```

because the `R`, `seeding`, and `infection_noise` component used the
defaults. No default is given for the generation time distribution,
since it is pathogen-specific.

💡 An advantage of coding each component explicitly is that your
modeling decisions are clearly documented.

### Estimation

We have assigned the five `EpiSewer` modules to individual variables,
which can now be passed to the `EpiSewer()` function to estimate the
effective reproduction number:

``` r
ww_measurements <- model_measurements()
options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
ww_result <- EpiSewer(
  data = sewer_data(measurements = measurements_sparse),
  measurements = ww_measurements,
  sampling = ww_sampling,
  sewage = ww_sewage,
  shedding = ww_shedding,
  infections = ww_infections,
  fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = 1000, iter_sampling = 1000, chains = 4), model = model_stan_opts(package = "EpiSewer")),
  run_fit = FALSE
)
```

Note that this time we don’t need a `data` or `assumptions` argument,
because we provided all data and assumptions to the respective module
components.

The above model is identical to the default model fitted in the
[README](../README.md), so we don’t repeat the results again.

Remember that you can print a summary of the full modeling details from
the job object:

``` r
ww_result$job$model
#> measurements
#>  |- concentrations_observe
#>  |- noise_estimate
#>  |- LOD_none
#> 
#> sampling
#>  |- sample_effects_none
#> 
#> sewage
#>  |- flows_observe
#>  |- residence_dist_assume
#> 
#> shedding
#>  |- incubation_dist_assume
#>  |- shedding_dist_assume
#>  |- load_per_case_assume
#>  |- load_variation_none
#> 
#> infections
#>  |- generation_dist_assume
#>  |- R_estimate_splines
#>  |- seeding_estimate_rw
#>  |- infection_noise_estimate
```

### Full model

For your convenience, we here provide the full model using the explicit
component specifications in one code block. This can serve as a starting
point to adapt the specification and try out different modeling options
(check out the function documentation for the module functions, or try
out the `component_functions()` helper).

Also, don’t forget that the modeling functions have further arguments to
customize priors and further modeling details. Again, the function
documentation is your friend.

``` r
ww_measurements <- model_measurements(
  concentrations = concentrations_observe(measurements = measurements_sparse),
  noise = noise_estimate(),
  LOD = LOD_none()
)

ww_sampling <- model_sampling(
  sample_effects = sample_effects_none()
)

ww_sewage <- model_sewage(
  flows = flows_observe(flows = data_zurich$flows),
  residence_dist = residence_dist_assume(residence_dist = c(1))
)

ww_shedding <- model_shedding(
  incubation_dist = incubation_dist_assume(get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4, maxX = 10)),
  shedding_dist = shedding_dist_assume(get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397, maxX = 30)),
  load_per_case = load_per_case_assume(6e+11),
  load_variation = load_variation_none()
)

ww_infections <- model_infections(
  generation_dist = generation_dist_assume(get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12)),
  R = R_estimate_splines(),
  seeding = seeding_estimate_rw(),
  infection_noise = infection_noise_estimate()
)

options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
ww_result <- EpiSewer(
  measurements = ww_measurements,
  sampling = ww_sampling,
  sewage = ww_sewage,
  shedding = ww_shedding,
  infections = ww_infections,
  fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = 1000, iter_sampling = 1000, chains = 4))
)
```
