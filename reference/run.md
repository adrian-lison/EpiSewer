# Fit an EpiSewer model.

This function allows to (re-)run the model fitting from an existing
EpiSewerJob object. It is especially useful when combined with the
`run_fit=FALSE` option in
[`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md).
This allows to first create the EpiSewerJob object and then run the
model fitting later or on a different machine.

## Usage

``` r
run(job, ...)
```

## Arguments

- job:

  An EpiSewerJob object as returned by
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md).
