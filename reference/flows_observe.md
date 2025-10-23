# Observe wastewater flows

This option accounts for daily wastewater flow volumes measured at the
sampling site. The flow can change due to rainfall or industrial
discharge, and directly influences pathogen concentrations in the
wastewater.

## Usage

``` r
flows_observe(
  flows = NULL,
  date_col = "date",
  flow_col = "flow",
  modeldata = modeldata_init()
)
```

## Arguments

- flows:

  A `data.frame` with each row representing one day. Must have at least
  a column with dates and a column with flow measurements.

- date_col:

  Name of the column containing the dates.

- flow_col:

  Name of the column containing the flows.

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

The flow volume unit should be the same as for the concentration
measurements, e.g. if concentrations are measured in gc/mL, then the
flow should be in mL as well.

## See also

Other flow models:
[`flows_assume()`](https://adrian-lison.github.io/EpiSewer/reference/flows_assume.md)
