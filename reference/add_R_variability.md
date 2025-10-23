# Add change point model for the variability of Rt

Add change point model for the variability of Rt

## Usage

``` r
add_R_variability(
  length_R,
  h,
  length_seeding,
  partial_window,
  partial_generation,
  changepoint_dist,
  modeldata
)
```

## Arguments

- length_R:

  Length of the modeled Rt time series

- changepoint_dist:

  How many days should the change points be apart? If zero, no change
  points are modeled (fixed R variability).

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

The change points are placed going backwards from the end of the time
series, i.e. on the last day, then `changepoint_dist` days before that,
and so on. The distance between the first and second change point at the
start of the time series (partly falls into the seeding phase) varies
between half and double the `changepoint_dist`.
