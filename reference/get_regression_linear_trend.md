# Obtain regression vector to estimate a linear trend

This function provides the solution vector from a (weighted) linear
regression with one independent variable x, where x may represent points
in a time series.

## Usage

``` r
get_regression_linear_trend(x, weights = rep(1/length(x), length(x)))
```

## Arguments

- x:

  A `vector` with observations of the single independent variable (for
  example points in time).

- weights:

  Optional, a `vector` with weights for each observation.

## Value

A `vector` which contains the relevant row of the regression solution
matrix that represents the trend coefficient. By multiplying this vector
with the observed time series values y, you directly get the (weighted
least squares) estimate of the trend.

## Examples

``` r
x = 1:6
y = x*4 + rnorm(6, 0, 0.1)
w <- c(0.1, 0.1, 0.1, 0.2, 0.2, 0.3) # weights
trend_reg <- get_regression_linear_trend(x, weights = w)
#> Error in get_regression_linear_trend(x, weights = w): could not find function "get_regression_linear_trend"
as.vector(trend_reg %*% y) # this gives you the trend estimate
#> Error: object 'trend_reg' not found
summary(lm(y ~ x, weights = w))$coefficients["x","Estimate"] # should be the same
#> [1] 3.986604
```
