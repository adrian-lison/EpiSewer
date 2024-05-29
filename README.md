
# EpiSewer: Estimate Reproduction Numbers from Wastewater Measurements <img src='man/figures/EpiSewer_logo.png' align="right" height="140" />

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![MIT
license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/epinowcast/epinowcast/blob/master/LICENSE.md/)
[![DOI](https://zenodo.org/badge/676114941.svg)](https://zenodo.org/doi/10.5281/zenodo.10569101)

The `EpiSewer` R package provides a Bayesian generative model to
estimate effective reproduction numbers from concentration measurements
at a wastewater treatment plant (or other sampling site) over time. This
allows to track the transmission dynamics of a pathogen in the
associated catchment population. The `EpiSewer` model is tailored to the
specifics of wastewater concentration measurements, offers comprehensive
uncertainty quantification via MCMC sampling in `stan`, and provides
easily configurable modeling components.

<img src="man/figures/weekly_zurich_animated.gif" width="100%" />

## Model highlights

**Measurements**  
⭐ Non-daily / missing measurements  
⭐ Multiple measurements (replicates) per sample  
⭐ Multi-day composite samples  
⭐ Accurate ddPCR noise model  
⭐ Limit of detection (LOD) model

**Sampling**  
⭐ Sample batch effects (e.g. weekday or age-of-sample effects)

**Sewer**  
⭐ Flow normalization  
⭐ Sewer residence time distributions

**Shedding**  
⭐ Incubation period and shedding load distributions  
⭐ Individual-level shedding load variation

**Infections**  
⭐ Stochastic infection model with overdispersion  
⭐ Flexible $R_t$ smoothing (random walk, exponential smoothing,
splines)

## Installing the package

The development version of `EpiSewer` can be installed from GitHub as
shown below. Please note that the package is still highly experimental
and may be subject to breaking changes.

``` r
remotes::install_github("adrian-lison/EpiSewer", dependencies = TRUE)
```

`EpiSewer` also requires CmdStan to be installed on your system. This
can be done using the `install_cmdstan()` function from `cmdstanr`. If
you experience any problems installing CmdStan, see the [cmdstanr
vignette](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) for help.

``` r
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan(cores = 2) # use more cores to speed up
```

## Introduction

This is a quick introduction to using the `EpiSewer` package. To learn
more about modeling with `EpiSewer`, see the [model
specification](vignettes/model-specification.md) and [detailed
example](vignettes/detailed-example.md) vignettes.

### Loading the package

`EpiSewer` comes with `data.table` for data manipulation and `ggplot2`
for plotting.

``` r
library(EpiSewer)
library(data.table)
library(ggplot2)
```

### Data

#### Load exemplary data

The below wastewater data are from Zurich, Switzerland. They are
provided by EAWAG (Swiss Federal Institute of Aquatic Science and
Technology) to the public domain (CC BY 4.0 license).

``` r
data_zurich <- SARS_CoV_2_Zurich
```

#### Inspect the data

The dataset contains `measurements`, `flows` and metadata about `units`.
There is also data about confirmed cases in the catchment area, which we
do not need explicitly, but can use to inform some of our prior
assumptions.

``` r
str(data_zurich)
#> List of 4
#>  $ measurements:Classes 'data.table' and 'data.frame':   120 obs. of  2 variables:
#>   ..$ date         : Date[1:120], format: "2022-01-01" "2022-01-02" ...
#>   ..$ concentration: num [1:120] NA NA 456 748 574 ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ flows       :Classes 'data.table' and 'data.frame':   120 obs. of  2 variables:
#>   ..$ date: Date[1:120], format: "2022-01-01" "2022-01-02" ...
#>   ..$ flow: num [1:120] 3.41e+11 3.41e+11 1.59e+11 1.61e+11 3.72e+11 ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ cases       :Classes 'data.table' and 'data.frame':   120 obs. of  2 variables:
#>   ..$ date : Date[1:120], format: "2022-01-01" "2022-01-02" ...
#>   ..$ cases: num [1:120] NA NA 1520 1727 1639 ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ units       :List of 2
#>   ..$ concentration: chr "gc/mL"
#>   ..$ flow         : chr "mL/day"
```

The measurements data contains daily viral concentration measurements
(in gc/mL, gc = gene copies) for the SARS-CoV-2 N1 gene at the
wastewater treatment plant in Zurich. Some days have missing
measurements, but this is no problem: `EpiSewer` naturally accounts for
missing values during estimation.

``` r
data_zurich$measurements
#>            date concentration
#>   1: 2022-01-01            NA
#>   2: 2022-01-02            NA
#>   3: 2022-01-03      455.7580
#>   4: 2022-01-04      747.8792
#>   5: 2022-01-05      573.7020
#>  ---                         
#> 116: 2022-04-26      182.3664
#> 117: 2022-04-27      379.2421
#> 118: 2022-04-28      439.7427
#> 119: 2022-04-29      394.0033
#> 120: 2022-04-30      293.1242
```

To show the handling of missing data more clearly, we make our data
artificially sparse by keeping only measurements that were made on
Mondays and Thursdays. This data only has two measurements per week, but
we can still use it to estimate reproduction numbers from it.

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

The flows data tracks the daily flow (in mL/day) at the treatment plant
in Zurich. The flow data will be used to normalize the concentration
measurements. This helps to account for environmental factors such as
rainfall. It is important that the flow uses the same volume unit as the
concentration (mL here in both cases).

``` r
data_zurich$flows
#>            date         flow
#>   1: 2022-01-01 341163082483
#>   2: 2022-01-02 341163082483
#>   3: 2022-01-03 158972454250
#>   4: 2022-01-04 160849294683
#>   5: 2022-01-05 372301227483
#>  ---                        
#> 116: 2022-04-26 330658518200
#> 117: 2022-04-27 206045734750
#> 118: 2022-04-28 158372697017
#> 119: 2022-04-29 152155673150
#> 120: 2022-04-30 208684746067
```

Note: In contrast to the concentration measurements, the flow data must
not have missing values. If flow information is missing for certain
days, make sure to impute it using a suitable method before passing it
to `EpiSewer`.

#### Gather the data

If the measurements data and flow data are stored in separate variables,
we can combine them using the helper function `sewer_data()`. We here
use our artificially sparse measurements.

``` r
ww_data <- sewer_data(measurements = measurements_sparse, flows = data_zurich$flows)
```

### Assumptions

In order to estimate the effective reproduction number from wastewater
concentration measurements, we must make a number of assumptions.

#### Distributions

- **Generation time distribution**: Distribution of the time between a
  primary infection and its resulting secondary infections
- **Incubation period distribution**: Time between infection and symptom
  onset. It is here assumed that shedding only starts with symptom
  onset.
- **Shedding load distribution**: Distribution of the load shed by an
  average individual over time since their symptom onset.

The generation time, incubation period and shedding load distribution
are all disease-specific and are typically obtained from literature.
`EpiSewer` requires these distributions to be discretized, and offers
specific functions to obtain discretized versions of popular continuous
probability distributions.

``` r
generation_dist <- get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12)
incubation_dist <- get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4, maxX = 10)
shedding_dist <- get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397, maxX = 30)
```

#### Shedding load per case

In addition to the shedding load distribution, we also require a scaling
factor that describes how much total viral material is shed by the
average infected individual and is detectable at the sampling site. We
call this scaling factor *shedding load per case* because it scales the
estimated number of cases/infections in our model. The shedding load per
case (gc/person) depends both on biological factors as well as on the
sewage system.

Fortunately, if we are only interested in the effective reproduction
number, a rough estimate for the shedding load per case is sufficient:
getting it wrong cannot strongly bias the reproduction number estimates,
and will distort the uncertainty of estimates only in rather extreme
cases.

If we also have case data available (ideally catchment-specific, or we
might down-scale case counts from a higher aggregation level to roughly
the size of our catchment), then we can use the function
`suggest_load_per_case` to obtain a rough estimate of a suitable
shedding load per case. The argument `ascertainment_prop` can be used to
account for underdetection of infections. For simplicity, we here assume
`ascertainment_prop=1`, meaning that 100% of infections become confirmed
cases (which is often not realistic).

``` r
suggest_load_per_case(
  data_zurich$measurements,
  data_zurich$cases,
  data_zurich$flows,
  ascertainment_prop = 1
)
#> [1] 1.3e+11
```

Hence we will assume `1e+11 gc/person` as the average shedding load.

``` r
load_per_case <- 1e+11
```

#### Combining all assumptions

The above assumptions can be conveniently recorded in a list by using
the function `sewer_assumptions()`.

``` r
ww_assumptions <- sewer_assumptions(
  generation_dist = generation_dist,
  incubation_dist = incubation_dist,
  shedding_dist = shedding_dist,
  load_per_case = load_per_case
)
```

### Estimation

Now that we have the data and the necessary assumptions, we can use
`EpiSewer` to estimate the effective reproduction number. We here use
the default model and settings provided by `EpiSewer`. With the helper
function `set_fit_opts()` we specify our sampling approach: we apply
Hamiltonian MCMC sampling via stan, using 4 chains with 500 warmup and
500 sampling iterations each.

Stan regularly provides updates about the progress of the sampler. The
overall runtime will depend on your hardware resources, the size of the
data, the complexity of the model used, and how well the model actually
fits the data.

``` r
options(mc.cores = 4) # allow stan to use 4 cores, i.e. one for each chain
ww_result <- EpiSewer(
  data = ww_data,
  assumptions = ww_assumptions,
  fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = 500, iter_sampling = 500, chains = 4))
)
```

💡 We are here using a shorthand version of the `EpiSewer()` function
which uses all available default settings. To find out how the
`EpiSewer` model can be customized, see the [model
specification](vignettes/model-specification.md) and [detailed
example](vignettes/detailed-example.md) vignettes.

### Plotting the results

Great, the sampling has completed! Now we can inspect the results of our
model. Plots are a convenient way to get a quick overview.

#### Model fit

It is good practice to first assess how well the model actually fitted
to the data. For this, we plot the observed concentration measurements
against the ones predicted by the model. The ribbon shows the 95%
credible interval. The black dots show measurements that we have
actually observed in our artificially sparse data set (only Mondays and
Thursdays). The grey crosses show all the other measurements that were
not available to the model. As we can see, the model still achieved a
decent fit on these missing days.

``` r
plot_concentration(ww_result, measurements = data_zurich$measurements)
```

<img src="man/figures/README-concentration-1.png" width="100%" /> By
default, the function `plot_concentration()` includes the modeled
measurement noise in the predicted concentrations. If you only want to
see the predicted expected concentration without noise, set
`include_noise = FALSE`.

``` r
plot_concentration(ww_result, measurements = data_zurich$measurements, include_noise = FALSE)
```

<img src="man/figures/README-concentration_noise-1.png" width="100%" />

#### Time-varying effective reproduction number

Since the model fit looked decent at first glance, we now inspect our
main parameter of interest, the effective reproduction number over time.
Again, the ribbon shows the 95% credible interval. We see that even with
only two measurements per week, we get a fairly clear signal for R being
above or below the critical threshold of 1.

``` r
plot_R(ww_result)
```

<img src="man/figures/README-R-1.png" width="100%" />

Note that the estimates for R go back further into the past than our
observations. This is due to the delay from infection to shedding,
i.e. concentration measurements observed today are mostly a signal of
infections in the past. This is also why the R estimates close to the
present are strongly uncertain. The measurements observed until the
present provide only a delayed signal about transmission dynamics, so
the most recent R estimates are informed by only little data.

#### Latent parameters

We can also inspect other internal parameters of the model. This may
give us a better understanding of the model fit and may help to diagnose
certain problems.

First, we plot the estimated load over time. This is the expected total
load that arrived at the sampling side on a given day. We can see this
is roughly on the order 1e14 to 3e14.

``` r
plot_load(ww_result)
```

<img src="man/figures/README-load-1.png" width="100%" /> Next, we plot
the estimated number of infections over time. The time series follows a
very similar trend as the load. Using our assumed shedding load per case
of 1e11, we can see that the peak load (roughly 3.3e14) corresponds to
3.3e14/1e11=3300 infections. The peak in infections shown below seems to
be a bit higher, and it is also a bit earlier. The reason for this
difference is that infected individuals only begin shedding after their
infection and then shed over a longer period of time (as defined by the
incubation period and shedding load distribution). This makes the load a
slightly delayed and blurred signal of the infections.

``` r
plot_infections(ww_result)
```

<img src="man/figures/README-infections-1.png" width="100%" /> We can
also compare the prior and posterior distribution for certain parameters
of the model to see how much they were informed by the data. For
example, we can inspect the coefficient of variation (CV) of the
observation noise:

``` r
plot_prior_posterior(ww_result, "measurement_noise_cv")
```

<img src="man/figures/README-prior_posterior_noise_cv-1.png" width="100%" />
As we can see, the prior (in grey) and posterior (in blue) differ
substantially, indicating that this parameter is well informed by the
data.

We can also inspect the CV for the individual-level variation in
shedding loads which is estimated by EpiSewer:

``` r
plot_prior_posterior(ww_result, "load_variation_cv")
```

<img src="man/figures/README-prior_posterior_load_cv-1.png" width="100%" />
As we can see here, the prior and posterior are almost identical, which
shows that this parameter cannot be really estimated from the data.
Thus, our uncertainty about individual-level shedding load variation is
not reduced by fitting to the data and solely informed by our broad
prior.

### More details

We can further inspect our results object. It has three attributes:

``` r
names(ww_result)
#> [1] "job"        "stan_model" "checksums"  "summary"    "fitted"
```

The `job` attribute stores all information about the job that was
defined via `EpiSewer`. It contains the data used, alongside
meta-information, and the settings for the sampler. By calling
`run(ww_result$job)`, the job can be run again.

``` r
names(ww_result$job)
#> [1] "job_name"      "jobarray_size" "data"          "model"        
#> [5] "init"          "fit_opts"      "priors_text"   "metainfo"     
#> [9] "overwrite"
```

In particular, we can print a concise summary of the modeling details
(the `EpiSewer` defaults in this case):

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
#>  |- load_variation_estimate
#> 
#> infections
#>  |- generation_dist_assume
#>  |- R_estimate_splines
#>  |- seeding_estimate_rw
#>  |- infection_noise_estimate (overdispersion = FALSE)
```

➡️ Check out the [model specification](vignettes/model-specification.md)
and [detailed example](vignettes/detailed-example.md) vignettes to
better understand the modeling details and customize the `EpiSewer`
model.

The `summary` attribute stores summarized results for important
parameters from the model.

``` r
names(ww_result$summary)
#> [1] "R"                           "R_samples"                  
#> [3] "expected_infections"         "expected_infections_samples"
#> [5] "infections"                  "infections_samples"         
#> [7] "expected_load"               "expected_concentration"     
#> [9] "concentration"
```

For example, we can access the exact estimates for the reproduction
number.

``` r
ww_result$summary$R
#>            date      mean    median lower_0.95 lower_0.5 upper_0.5 upper_0.95
#>   1: 2021-12-06 1.0677532 1.0597600  0.7376447 0.9644702  1.168880   1.392586
#>   2: 2021-12-07 1.0686402 1.0612000  0.7641206 0.9684573  1.163675   1.383842
#>   3: 2021-12-08 1.0696771 1.0627000  0.7790942 0.9695755  1.162902   1.379000
#>   4: 2021-12-09 1.0708548 1.0631000  0.7895383 0.9734535  1.161805   1.367851
#>   5: 2021-12-10 1.0721631 1.0644950  0.8011415 0.9788717  1.161825   1.358601
#>  ---                                                                         
#> 140: 2022-04-24 0.9772645 0.9715175  0.7512149 0.8977970  1.052825   1.226163
#> 141: 2022-04-25 0.9771118 0.9699495  0.7359722 0.8900868  1.056298   1.250144
#> 142: 2022-04-26 0.9771515 0.9729495  0.7167654 0.8845507  1.061923   1.281037
#> 143: 2022-04-27 0.9774270 0.9738840  0.6972232 0.8773860  1.070855   1.309254
#> 144: 2022-04-28 0.9778231 0.9760445  0.6676156 0.8707375  1.075880   1.340137
#>      seeding
#>   1:    TRUE
#>   2:    TRUE
#>   3:    TRUE
#>   4:    TRUE
#>   5:    TRUE
#>  ---        
#> 140:   FALSE
#> 141:   FALSE
#> 142:   FALSE
#> 143:   FALSE
#> 144:   FALSE
```

The `fitted` attribute provides access to all details of the fitted stan
model. See
[cmdstanr](https://mc-stan.org/cmdstanr/reference/CmdStanMCMC.html) for
details.

``` r
names(ww_result$fitted)
#>  [1] ".__enclos_env__"            "functions"                 
#>  [3] "runset"                     "num_chains"                
#>  [5] "inv_metric"                 "diagnostic_summary"        
#>  [7] "sampler_diagnostics"        "loo"                       
#>  [9] "clone"                      "draws"                     
#> [11] "output"                     "initialize"                
#> [13] "code"                       "profiles"                  
#> [15] "return_codes"               "metadata"                  
#> [17] "time"                       "data_file"                 
#> [19] "latent_dynamics_files"      "profile_files"             
#> [21] "output_files"               "save_data_file"            
#> [23] "save_profile_files"         "save_latent_dynamics_files"
#> [25] "save_output_files"          "cmdstan_diagnose"          
#> [27] "cmdstan_summary"            "summary"                   
#> [29] "lp"                         "constrain_variables"       
#> [31] "variable_skeleton"          "unconstrain_draws"         
#> [33] "unconstrain_variables"      "hessian"                   
#> [35] "grad_log_prob"              "log_prob"                  
#> [37] "init_model_methods"         "init"                      
#> [39] "save_object"                "expose_functions"          
#> [41] "print"                      "num_procs"
```

For example, we can use this to show sampler diagnostics for each chain:

``` r
ww_result$fitted$diagnostic_summary()
#> $num_divergent
#> [1] 0 0 0 0
#> 
#> $num_max_treedepth
#> [1] 0 0 0 0
#> 
#> $ebfmi
#> [1] 0.8404351 0.8599965 0.7253768 0.8287051
```

Finally, the `checksums` attribute gives us several checksums that
uniquely identify the job which was run. These can be used to check
whether two job results used the same or different models, input data,
fitting options, or inits. If all checksums are identical (and the seed
is not NULL), then the results should also be identical.

``` r
ww_result$checksums
#> $model
#> [1] "77c97cfa2bfba2be740b2fbae8ae59b1"
#> 
#> $input
#> [1] "079bf4474da23e47636126383424f64c"
#> 
#> $fit_opts
#> [1] "3d0f79a45e43b33e7c77ffaee9760d1c"
#> 
#> $init
#> [1] "4c66acdaea245b71fe10c9366a583479"
```
