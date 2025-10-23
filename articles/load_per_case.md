# The "shedding load per case" assumption

In addition to the shedding load distribution, we also require a scaling
factor that describes how much total genetic material is shed by the
average infected individual and is detectable at the sampling site. We
call this scaling factor *shedding load per case* because it scales the
estimated number of infections in our model. The shedding load per case
(gc/person) depends both on biological factors as well as on the sewage
system and quantification method.

Fortunately, if we are only interested in the effective reproduction
number, a rough estimate for the shedding load per case is sufficient:
getting it wrong cannot strongly bias the reproduction number estimates,
and will distort the uncertainty of estimates only in rather extreme
cases.

## Disclaimer on interpretation of estimated infections

One important aspect to keep in mind when interpreting the number of
infections modeled in EpiSewer is that when modeling wastewater data,
the reference population for the infection dynamics is often not clear:

- The population in the catchment does not have to be identical with the
  population relevant for the infection process. For example, if the
  catchment area of the sampled treatment plant covers only a certain
  region of a city, then the population relevant for shedding might only
  be a subset of the overall population in which relevant transmissions
  occur.
- The population in the catchment may also vary over time due to
  every-day mobility or special events.
- The same probably applies even to catchments covering a whole city if
  there is sufficient mixing with populations from other cities.

Ultimately, this means that the population modeled in the infection
process might be substantially larger than the catchment population.
EpiSewer currently has no specific functionality to represent this
difference (although it is planned for future versions).

This is also another reason why the number of infections should not be
interpreted in terms of true absolute incidence or prevalence (in
addition to the fact that it ultimately depends on the assumed load per
case, which is very difficult to specify correctly).

## Calibration to minimum number of infections

By default, `EpiSewer` will calibrate the load per case such that the
estimated number of infections does not go considerably below 10
infections per day. This is a robust and conservative default in the
sense that it most likely *underestimates* the true number of infections
(see considerations above) - which in turns means that our uncertainty
of epidemiological parameters is potentially too large, but not too
small. You can adapt this default with the `min_cases` argument in
[`sewer_assumptions()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_assumptions.md):

``` r
ww_assumptions <- sewer_assumptions(
  min_cases = 30
  # add further assumptions here
)
```

or using the
[`load_per_case_calibrate()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_calibrate.md)
function:

``` r
load_per_case_calibrate(min_cases = 30)
```

    ## shedding
    ##  |- load_per_case_calibrate

❗ Note that the default is deliberately not `min_cases=1`. This is
because the continuous approximations used by the stochastic infection
model of EpiSewer start to break down when infections numbers become
extremely small. This can affect sampling speed and accuracy of
estimates. Aside from this, it is also rather unlikely for most
catchments to detect genetic material of a pathogen if it is only shed
by a single individual in the whole catchment.

## Calibration to case data

If we also have case data available, we can provide this to EpiSewer via

``` r
ww_data <- sewer_data(
  cases = SARS_CoV_2_Zurich$cases
  # add further data here
  )
```

or using the
[`load_per_case_calibrate()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_calibrate.md)
function:

``` r
load_per_case_calibrate(cases = SARS_CoV_2_Zurich$cases)
```

    ## shedding
    ##  |- load_per_case_calibrate

This will calibrate the shedding load per case such that the estimated
number of infections roughly matches the observed number of cases.

If you want to do and check this manually, you can use the function
[`suggest_load_per_case()`](https://adrian-lison.github.io/EpiSewer/reference/suggest_load_per_case.md).
The argument `ascertainment_prop` can be used to account for
underdetection of infections. For simplicity, we here assume
`ascertainment_prop=1`, meaning that 100% of infections become confirmed
cases (which is often not realistic).

``` r
suggest_load_per_case(
  SARS_CoV_2_Zurich$measurements,
  SARS_CoV_2_Zurich$cases,
  SARS_CoV_2_Zurich$flows,
  ascertainment_prop = 1
)
```

    ## [1] 1.3e+11
