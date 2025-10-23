# Dampen forecasts

This option dampens the forecast of Rt so that trends in transmission
will level off after some time. This prevents unrealistic extrapolation
of transmission dynamics.

## Usage

``` r
damping_assume(damping = 0.95, modeldata = modeldata_init())
```

## Arguments

- damping:

  The forecast damping parameter. A value of 1 means no damping, a value
  of 0 means flat forecast. The default is 0.8, which levels off after a
  horizon of approximately 2 weeks.

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

The applied damping is exponential, i.e. the trend is reduced by
`damping^1` on the first forecast day, by `damping^2` on the second
forecast day, and so on.
