# Solve Dirichlet Parameters for Soft Changepoint Model

This function solves for the parameters p_i of a Dirichlet distribution
given a total sum alpha and number of parameters n.

## Usage

``` r
changepoint_dirichlet_prior(
  alpha,
  n,
  min_break_dist = 1,
  prior_previous = NULL
)
```

## Arguments

- alpha:

  Total sum of all p_i

- n:

  Number of parameters to solve for

- min_break_dist:

  Minimum distance between changepoints. For first segment, can be left
  empty/set to 1.

- prior_previous:

  Vector with alpha values of prior for previous segment. Required if
  prior_previous\>1.

## Value

A vector of p_i values
