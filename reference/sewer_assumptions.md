# Specify modeling assumptions

Specify model assumptions such as the generation time distribution or
shedding load distribution. This is a convenience function to collect
all assumptions in one object.

## Usage

``` r
sewer_assumptions(
  generation_dist = NULL,
  shedding_dist = NULL,
  shedding_reference = NULL,
  incubation_dist = NULL,
  min_cases = 10,
  load_per_case = NULL,
  residence_dist = c(1),
  ...
)
```

## Arguments

- generation_dist:

  Generation time distribution. The intrinsic distribution of the time
  between infection of a primary case and infection of its secondary
  cases. Will be automatically passed to
  [`generation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/generation_dist_assume.md).

- shedding_dist:

  Shedding load distribution. Describes how the total load shed by an
  individual is distributed over time (and therefore sums to 1). Will be
  automatically passed to
  [`shedding_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/shedding_dist_assume.md).

- shedding_reference:

  Is the shedding load distribution relative to the day of `"infection"`
  or the day of `"symptom_onset"`? This is important because shedding
  load distributions provided in the literature are sometimes by days
  since infection and sometimes by days since symptom onset. If
  `shedding_reference="symptom_onset"`, EpiSewer also needs information
  about the incubation period distribution. Will be automatically passed
  to
  [`shedding_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/shedding_dist_assume.md).

- incubation_dist:

  Incubation period distribution. The incubation period is the time
  between infection and symptom onset. This assumption is used when
  `shedding_reference="symptom_onset"`, i.e. to support shedding load
  distributions referenced by days since symptom onset. Will be
  automatically passed to
  [`incubation_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/incubation_dist_assume.md).

- min_cases:

  The minimum incidence (new cases per day) during the modeled time
  period (default is 10). This assumption is used to calibrate the
  `load_per_case` factor (a scaling factor that describes how many
  pathogen particles are shed by the average infected individual overall
  and how much of this is detectable at the sampling site). If
  `min_cases` is specified, `load_per_case` is chosen such that the
  estimated time series of infections does not go significantly below
  `min_cases`. Will be automatically passed to
  [`load_per_case_calibrate()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_calibrate.md).

- load_per_case:

  This argument allows to directly specify the average total load per
  case, instead of calibrating it to case data. Note that the
  `load_per_case` depends both on biological factors as well as on the
  specific sewage system and laboratory quantification. Important: To
  actually use this assumption, the model option
  [`load_per_case_assume()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_assume.md)
  must be specified in
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  via
  [`model_shedding()`](https://adrian-lison.github.io/EpiSewer/reference/model_shedding.md).
  If you want to use this option to do a manual calibration to case
  data, the helper function
  [`suggest_load_per_case()`](https://adrian-lison.github.io/EpiSewer/reference/suggest_load_per_case.md)
  might also be useful.

- residence_dist:

  Sewer residence time distribution for pathogen particles. By default,
  `EpiSewer` assumes that particles arrive at the sampling site within
  the day of shedding. However, for larger sewage systems, particles may
  travel longer than a day depending on where and when they were shed
  into the wastewater. Will be automatically passed to
  [`residence_dist_assume()`](https://adrian-lison.github.io/EpiSewer/reference/residence_dist_assume.md).

- ...:

  Further assumptions to be passed to
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md).

## Value

A `list` with all assumptions supplied. Can be passed to the
`assumptions` argument in
[`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md).

## Details

Note that if case data is supplied via `sewer_data`, the `min_cases`
assumption is overwritten and the supplied case data is used for a more
accurate calibration instead.
