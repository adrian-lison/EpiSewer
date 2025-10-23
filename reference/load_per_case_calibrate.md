# Calibrate the average load per case using case count data

This option calibrates the average total shedding load per case based on
the relationship between measured concentrations and case counts. The
goal is to obtain a `load_per_case` that is on the right order of
magnitude - this will not be sufficient for accurate prevalence
estimation from wastewater, but fully sufficient for monitoring trends
and estimating Rt.

## Usage

``` r
load_per_case_calibrate(
  cases = NULL,
  min_cases = NULL,
  ascertainment_prop = 1,
  measurement_shift = seq(-7, 7),
  shift_weights = 1/(abs(measurement_shift) + 1),
  date_col = "date",
  case_col = "cases",
  signif_fig = 2,
  modeldata = modeldata_init()
)
```

## Arguments

- cases:

  A `data.frame` of case numbers with each row representing one day.
  Must have at least a column with dates and a column with case numbers.

- min_cases:

  This is an alternative to supplying a `data.frame` of cases. If
  `min_cases` is specified, then `load_per_case` is calibrated based on
  the smallest detected load in the wastewater. EpiSewer uses
  `min_cases = 10` as a default if no case data is supplied (see
  [`sewer_assumptions()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_assumptions.md)).

- ascertainment_prop:

  Proportion of all cases that get detected / reported. Can be used to
  account for underreporting of infections. Default is
  `ascertainment_prop=1`, meaning that 100% of infections become
  confirmed cases.

- measurement_shift:

  The specific timing between wastewater concentrations and case numbers
  depends on reporting delays and shedding profiles and is typically
  uncertain. This argument allows to shift the concentration and case
  number time series relative to each other and to average over several
  potential lags/leads, as specified by an integer vector. The default
  is `measurement_shift = seq(-7,7)`, i.e. a shift of concentrations
  between up to one week before and after case numbers.

- shift_weights:

  Weights for the shifted comparisons. Must be an numeric vector of the
  same length as `measurement_shift`. If `NULL` (default), the weights
  are chosen to be approximately inversely proportional to the shift
  distance.

- date_col:

  Name of the date column in the provided cases `data.frame`.

- case_col:

  Name of the column containing the case numbers.

- signif_fig:

  Significant figures to round to. Since this heuristic only provides
  crude estimates which should not be over-interpreted, the result gets
  rounded (this also increases stability when fitting the model at
  different points in time). Default is rounding to the 2 most
  significant figures.

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

In `EpiSewer`, the `load_per_case` serves as a scaling factor describing
how many pathogen particles are shed by the average infected individual
overall and how much of this is detectable at the sampling site. This
depends both on biological factors as well as on the specific sewage
system and laboratory quantification. It is therefore almost always
necessary to assume the load per case based on a comparison of measured
concentrations/loads and case numbers.

If a `data.frame` of cases is supplied via the `cases` argument, the
average load per case is determined by fitting a linear regression model
with loads (computed using concentrations and flows) as dependent
variable and case numbers as independent variable. This does not
explicitly account for shedding profiles or reporting delays, but the
`measurement_shift` argument allows to average over a set of relative
shifts between the two time series.

## See also

Other load per case functions:
[`load_per_case_assume()`](https://adrian-lison.github.io/EpiSewer/reference/load_per_case_assume.md)
