# Define a normal prior on the log scale in modeldata using median and factor

This parameterization is useful to define a best guess for the value of
a parameter (median), and a maximum factor by which we expect to deviate
from that guess.

## Usage

``` r
set_prior_normal_log(
  param,
  unit_median = NULL,
  unit_q5 = NULL,
  unit_q95 = NULL,
  unit_factor = NULL,
  paste_dist = ""
)
```

## Arguments

- param:

  Name of the parameter for which the prior is defined.

- unit_median:

  Median on the unit/natural scale.

- unit_q5:

  5% quantile of the distribution on the unit scale.

- unit_q95:

  95% quantile of the distribution on the unit scale.

- unit_factor:

  By which factor do we expect the true parameter value to differ at
  most from our prior median? For example, `unit_factor = 2` would mean
  that we expect the true parameter to be at most twice our prior
  median, and at least half of our prior median. This uses the two-sigma
  rule-of-thumb. Must be specified together with `unit_median.`

- paste_dist:

  Additional information that should be pasted to the dist description.

## Value

Prior specification for modeldata.

## Details

A normal prior on the log scale is effectively a log-normal prior on the
unit/natural scale.
