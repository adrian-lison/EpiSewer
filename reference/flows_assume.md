# Assume a constant wastewater flow

This option assumes a constant flow of wastewater at the sampling site.
Can be used as an approximation if no regular flow measurements are
available.

## Usage

``` r
flows_assume(flow_constant, modeldata = modeldata_init())
```

## Arguments

- flow_constant:

  Daily wastewater flow volume, assumed to be the same for all days.

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

Note that when the flow is unknown and therefore some arbitrary value
(e.g. `flow_constant=1`) is assumed, this must also be accounted for in
the assumed `load_per_case`. For this purpose, the function
[`suggest_load_per_case()`](https://adrian-lison.github.io/EpiSewer/reference/suggest_load_per_case.md)
offers a `flow_constant` argument.

## See also

Other flow models:
[`flows_observe()`](https://adrian-lison.github.io/EpiSewer/reference/flows_observe.md)
