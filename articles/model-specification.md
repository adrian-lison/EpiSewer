# Modeling in EpiSewer

This is a conceptual overview over the modeling functionality in
`EpiSewer`. See
[`vignette("model-definition")`](https://adrian-lison.github.io/EpiSewer/articles/model-definition.md)
for a mathematical definition of the underlying generative model, and
[`vignette("detailed-example")`](https://adrian-lison.github.io/EpiSewer/articles/detailed-example.md)
for an example vignette.

## Modules

`EpiSewer` uses 5 different modules to describe the data generating
process behind the wastewater measurements: `infections`, `shedding`,
`sewage`, `sampling`, and `measurements`. There is a 6th module to
specify `forecast` functionality. Each of these modules consists of a
number of module components, as shown below.

![The 6 modules of the EpiSewer
model](figures/specification-model-1.png)

The modules are defined using their corresponding module function,
i.e. by calling
[`model_infections()`](https://adrian-lison.github.io/EpiSewer/reference/model_infections.md),
[`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md),
[`model_sewage()`](https://adrian-lison.github.io/EpiSewer/reference/model_sewage.md),
[`model_sampling()`](https://adrian-lison.github.io/EpiSewer/reference/model_sampling.md),
[`model_measurements()`](https://adrian-lison.github.io/EpiSewer/reference/model_measurements.md),
or
[`model_forecast()`](https://adrian-lison.github.io/EpiSewer/reference/model_forecast.md).

## Modeling functions

Components in a module can be specified using suitable modeling
functions. There are 5 types of modeling functions:

- `_observe`: We provide observation data for this component. For
  example, we can use
  [`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md)
  if we have observed concentration measurements and want to fit the
  model to them.
- `_assume`: We assume the values for this component. For example, we
  can use
  [`generation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/generation_dist_assume.md)
  to provide a generation time distribution from the literature.
- `_calibrate`: This is similar to `_assume`, but instead of directly
  specifying the value for an assumption, we calibrate it to some other
  assumption or data. For example, we can use
  [`load_per_case_calibrate()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_calibrate.md)
  to calibrate the shedding load per case to case data (so that the
  estimated infections will roughly match the observed case numbers).
- `_estimate`: We estimate this component as a parameter of the model.
  For example, we can use
  [`noise_estimate()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate.md)
  if we don’t know how much noise the measurements have and want to
  estimate this from the data.
- `_none`: We do not model this component. For example, we can use
  [`sample_effects_none()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_none.md)
  if we don’t want to model any sample effects on the measurements.

To find out which modeling functions are available for a given
component, you can consult the documentation or use the helper
[`component_functions()`](https://adrian-lison.github.io/EpiSewer/reference/component_functions.md):

``` r
EpiSewer::component_functions("infection_noise")
#> [1] "infection_noise_none()"     "infection_noise_estimate()"
```

#### 💡 Multiple modeling options

Some components have multiple versions of the same modeling function
type. For example, there are several approaches to estimate the
reproduction number, such as
[`R_estimate_gp()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_gp.md)
(Gaussian process),
[`R_estimate_rw()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_rw.md)
(random walk),
[`R_estimate_ets()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_ets.md)
(exponential smoothing),
[`R_estimate_splines()`](https://adrian-lison.github.io/EpiSewer/reference/R_estimate_splines.md)
(smoothing splines) and more.

``` r
EpiSewer::component_functions("R")
#> [1] "R_estimate_gp()"                  "R_estimate_ets()"                
#> [3] "R_estimate_smooth_derivative()"   "R_estimate_piecewise()"          
#> [5] "R_estimate_changepoint_splines()" "R_estimate_approx()"             
#> [7] "R_estimate_rw()"                  "R_estimate_splines()"
```

#### ❗ Modeling restrictions

Not all components support all modeling types. For example, `EpiSewer`
currently only offers
[`generation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/generation_dist_assume.md),
but not `generation_dist_estimate`, `generation_dist_calibrate`, or
`generation_dist_observe`.

``` r
EpiSewer::component_functions("generation_dist")
#> [1] "generation_dist_assume()"
```

This is because estimating the generation time distribution from data is
not yet supported (but may be added in the future).

## Data and assumptions

The
[`sewer_data()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_data.md)
and
[`sewer_assumptions()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_assumptions.md)
functions are convenience functions that allow you to collect all
observation data and modeling assumptions in one place. This can improve
overview and allows to run `EpiSewer` under its default settings without
repeating the component definitions:

``` r
ww_data <- sewer_data(
  measurements = SARS_CoV_2_Zurich$measurements,
  flows = SARS_CoV_2_Zurich$flows,
  cases = SARS_CoV_2_Zurich$cases # optional
  )

ww_assumptions <- sewer_assumptions(
  generation_dist = get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4),
  shedding_dist = get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397),
  shedding_reference = "symptom_onset",
  incubation_dist = get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4),
)

EpiSewer(
  data = ww_data,
  assumptions = ww_assumptions
)
```

What happens under the hood is that when individual model components are
not provided with the data or assumptions they need, they search the
data and assumption arguments passed to the
[`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
function.

💡 It is always possible to mix both approaches and specify data or
assumptions explicitly in the model component:

``` r
# Leave out flows from the data
ww_data <- sewer_data(
  measurements = SARS_CoV_2_Zurich$measurements,
  #flows = SARS_CoV_2_Zurich$flows
  cases = SARS_CoV_2_Zurich$cases
)

# Leave out the generation time distribution
ww_assumptions <- sewer_assumptions(
  #generation_dist = get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4),
  shedding_dist = get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397),
  shedding_reference = "symptom_onset",
  incubation_dist = get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4),
)

# Provide flows directly to sewage module
ww_sewage <- model_sewage(
  flows = flows_observe(SARS_CoV_2_Zurich$flows)
)

# Provide generation time distribution directly to infections module
ww_infections <- model_infections(
  generation_dist = generation_dist_assume(
    get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4)
  )
)

# Combine everything
result <- EpiSewer(
  data = ww_data,
  assumptions = ww_assumptions,
  sewage = ww_sewage,
  infections = ww_infections
)
```

Note that if the same data or assumptions are supplied via the
[`sewer_data()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_data.md)/[`sewer_assumptions()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_assumptions.md)
*and* the individual component,
[`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
will compare both arguments and throw an error if they differ.
