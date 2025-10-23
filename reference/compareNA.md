# Compare entries of two vectors with potentially missing values

Compare entries of two vectors with potentially missing values

## Usage

``` r
compareNA(v1, v2)
```

## Arguments

- v1:

  First vector

- v2:

  Second vector, should be of same length as first one

## Value

A boolean vector with element-wise comparisons. Comparisons including
NAs are coded as FALSE.
