# Compute sigma parameter of Log-Normal from other quantities.

Compute sigma parameter of Log-Normal from other quantities.

## Usage

``` r
get_lognormal_sigma_alternative(
  mu,
  unit_mean = NULL,
  unit_sd = NULL,
  unit_q5 = NULL,
  unit_q95 = NULL
)
```

## Arguments

- mu:

  Mu parameter of Log-Normal distribution.

- unit_q5:

  5% quantile of distribution.

- unit_q95:

  95% quantile of distribution.

## Value

Sigma parameter of Log-Normal distribution.

## Details

Currently, only conversion from mu + a quantile is supported, but other
alternatives may be added.
