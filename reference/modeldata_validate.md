# Ensure that modeldata contains all necessary variables

Throws an error if necessary variables are missing, and replaces
optional variables with defaults indicating their absence in the stan
model.

## Usage

``` r
modeldata_validate(
  modeldata,
  data = list(),
  assumptions = list(),
  defaults = modeldata_defaults()
)
```

## Arguments

- modeldata:

  A `list` with all data variables (including priors) to be passed on to
  stan, alongside additional meta information and descriptions.

- data:

  A list with observation data as provided by
  [`sewer_data()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_data.md).

- assumptions:

  A list with assumptions as provided by
  [`sewer_assumptions()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_assumptions.md).

- defaults:

  A `list` with default values to be used for modeldata variables if not
  supplied in modeldata. For example, `numeric(0)` will often be
  supplied for optional parameters.

## Value

If no error is thrown due to missing mandatory variables, the same
modeldata object is returned again, where optional variables have been
replaced by zero-length dimension defaults for stan.

## Details

The reason why defaults are inserted for missing optional variables is
that stan does not explicitly allow default variables, so they are
specified with a little hack by using zero-length dimensions.
