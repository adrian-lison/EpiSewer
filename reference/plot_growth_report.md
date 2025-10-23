# Plot a growth report

For a given date, the growth report shows the probability that
infections have been growing constantly for at least 3, 7, 14, 21, and
28 days, respectively. The report uses a diverging bar plot which is
scaled between "very unlikely" (0% posterior probability) and "very
likely" (100% posterior probability).

## Usage

``` r
plot_growth_report(result, date = NULL, partial_prob = 0.8)
```

## Arguments

- result:

  Results object returned by
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  after model fitting. In contrast to other plotting functions, this
  cannot be a list of multiple result objects, because growth report
  plots are always for a single model fit.

- date:

  The date for which the growth report should be plotted. If NULL, the
  most recent date for which a reliable report can be provided is
  automatically selected, see `partial_prob`. Note that the date is
  backward looking. For example, if `date="2022-03-01"`, you will get
  the probability that have been growing in the last 3, 7, 14, 21, and
  28 days before this date.

- partial_prob:

  To select the most recent reliable date, we subtract a certain
  quantile of the shedding load distribution from the current date. For
  example, if `partial_prob=0.8` (default), we select the date for which
  80% of the shedding load of individuals infected before this date has
  been shed.

## Value

A growth report plot showing the probability that infections have been
growing without interruption for at least 3, 7, 14, 21, and 28 days,
respectively. Can be further manipulated using `ggplot2` functions to
adjust themes and scales, and to add further geoms.
