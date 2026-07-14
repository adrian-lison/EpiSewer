# Run a quick model fitting test with minimal iterations.

Modifies sampler settings to use very few iterations and then calls
[`run()`](https://adrian-lison.github.io/EpiSewer/reference/run.md).
Useful for quickly checking that a model specification is valid and that
the model can be fitted without errors.

## Usage

``` r
test_run(job)
```

## Arguments

- job:

  An `EpiSewerJob` or `EpiSewerJobResult` object.

## Value

An `EpiSewerJobResult` object.
