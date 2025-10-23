# Observe concentration measurements

This option fits the `EpiSewer` model to pathogen concentrations
measured in wastewater samples. It is suitable for different
quantification methods such as qPCR or dPCR. The measured concentrations
are by default modeled via a gamma likelihood.

## Usage

``` r
concentrations_observe(
  measurements = NULL,
  composite_window = 1,
  distribution = "gamma",
  date_col = "date",
  concentration_col = "concentration",
  replicate_col = NULL,
  n_averaged = 1,
  n_averaged_col = NULL,
  total_partitions_col = NULL,
  modeldata = modeldata_init()
)
```

## Arguments

- measurements:

  A `data.frame` with each row representing one measurement. Must have
  at least a column with dates and a column with concentration
  measurements.

- composite_window:

  Over how many days has each measured sample been collected? If 1
  (default), samples represent single days. If larger than 1, samples
  are assumed to be equivolumetric composites over several dates. In
  this case, the supplied dates represent the last day included in each
  sample.

- distribution:

  Parametric distribution for concentration measurements. Currently
  supported are "gamma" (default and recommended), "log-normal",
  "truncated normal", and "normal". The "truncated normal" and "normal"
  options are not recommended for use in practice.

- date_col:

  Name of the column containing the dates.

- concentration_col:

  Name of the column containing the measured concentrations.

- replicate_col:

  Name of the column containing the replicate ID of each measurement.
  This is used to identify multiple measurements made of a sample from
  the same date. Should be `NULL` if only one measurement per date was
  made.

- n_averaged:

  The number of replicates over which the measurements have been
  averaged. This is typically used as an alternative to providing
  several replicates per sample (i.e. the concentration provided in the
  `measurements` `data.frame` is the average of several replicates). Can
  be either a single number (it is then assumed that the number of
  averaged replicates is the same for each observation) or a vector (one
  value for each observation).

- n_averaged_col:

  Name of the column in the `measurements` data.frame containing the
  number of replicates over which the measurements have been averaged.
  This is an alternative to specifying `n_averaged`.

- total_partitions_col:

  Name of the column in the `measurements` data.frame containing the
  total number of partitions (e.g. droplets for ddPCR) in the dPCR
  reaction of each measurement. Only applies to concentration
  measurements obtain via dPCR. Can be used by the
  [`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
  and
  [`LOD_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
  modeling components. Note that this is really the *total* number of
  partitions, not just the number of positive partitions.

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
[`partitions_observe()`](https://adrian-lison.github.io/EpiSewer/reference/partitions_observe.md)
