# Observe positive dPCR partition counts

This option fits the `EpiSewer` model to positive partition counts in
digital PCR (dPCR), e.g. positive droplets in ddPCR, for the pathogen
target of interest. This allows the use of a dPCR-specific likelihood
using a Poisson distribution. For a more generic likelihood, see
[`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md).

## Usage

``` r
partitions_observe(
  data,
  composite_window = 1,
  date_col = "date",
  positive_partitions_col = "positive_partitions",
  total_partitions_col = "total_partitions",
  replicate_col = NULL,
  modeldata = modeldata_init()
)
```

## Arguments

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

## See also

Other observation types:
[`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md)
