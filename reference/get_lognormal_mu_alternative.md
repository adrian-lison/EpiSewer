# Compute mu parameter of Log-Normal from other quantities.

Compute mu parameter of Log-Normal from other quantities.

## Usage

``` r
get_lognormal_mu_alternative(
  unit_mean = NULL,
  unit_sd = NULL,
  unit_q5 = NULL,
  unit_q95 = NULL
)
```

## Arguments

- unit_q5:

  5% quantile of distribution.

- unit_q95:

  95% quantile of distribution.

## Value

Mu parameter of Log-Normal distribution.

## Details

Currently, only conversion from quantiles is supported, but other
alternatives may be added.
