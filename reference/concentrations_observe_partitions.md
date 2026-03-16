# Observe positive dPCR partition counts

This option fits the `EpiSewer` model to positive partition counts in
digital PCR (dPCR), e.g. positive droplets in ddPCR. This allows the use
of a dPCR-specific likelihood using a Binomial model for the number of
positive partitions observed. For a more generic likelihood, see
[`concentrations_observe()`](https://adrian-lison.github.io/EpiSewer/reference/concentrations_observe.md).

## Usage

``` r
concentrations_observe_partitions(
  measurements = NULL,
  composite_window = 1,
  date_col = "date",
  concentration_col = "concentration",
  positive_partitions_col = "positive_partitions",
  replicate_col = NULL,
  n_averaged = 1,
  n_averaged_col = NULL,
  total_partitions_col,
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

- date_col:

  Name of the column containing the dates.

- concentration_col:

  Name of the column containing the measured concentrations.

- positive_partitions_col:

  Name of the column in the `measurements` data.frame containing the
  number of positive partitions (e.g. positive droplets for ddPCR) in
  the dPCR reaction of each measurement. If several technical replicates
  are used, this should be the AVERAGE number of positive partitions per
  replicate.

- replicate_col:

  Name of the column containing the replicate ID of each measurement.
  This is used to identify multiple measurements made of a sample from
  the same date. Should be `NULL` if only one measurement per date was
  made.

- n_averaged:

  The number of technical replicates (i.e. repeated PCR runs) used for
  each sample or biological replicate. The concentration provided in the
  `measurements` `data.frame` is assumed to be the average or pooled
  estimate from several technical replicates. Can be either a single
  number (it is then assumed that the number of averaged replicates is
  the same for each observation) or a vector (one value for each
  observation).

- n_averaged_col:

  Name of the column in the `measurements` data.frame containing the
  number of technical replicates over which the measurements have been
  averaged/pooled. This is an alternative to specifying `n_averaged`.

- total_partitions_col:

  Name of the column in the `measurements` data.frame containing the
  number of total partitions (e.g. droplets for ddPCR) in the dPCR
  reaction of each measurement. If several technical replicates are
  used, this should be the AVERAGE number of valid partitions per
  replicate. Only applies when modeling concentration measurements via
  the dPCR-specific noise model. Can be used by the
  [`noise_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/noise_estimate_dPCR.md)
  and
  [`LOD_estimate_dPCR()`](https://adrian-lison.github.io/EpiSewer/reference/LOD_estimate_dPCR.md)
  modeling components. Note that this is really the number of *valid*
  partitions, not the number of positive partitions.

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
