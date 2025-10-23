# Lazy-execute a function once all necessary inputs are available

Lazy-execute a function once all necessary inputs are available

## Usage

``` r
tbc(
  f_name,
  f_expr,
  required = c(),
  required_values = NULL,
  advice = NULL,
  modeldata = NULL,
  calling_env = rlang::caller_env()
)
```

## Arguments

- f_name:

  Name of the "functions" that will be executed once all necessary
  inputs are available.

- f_expr:

  Expression with arbitrary R code in which attributes are assigned to
  `modeldata.`

- required:

  A character vector with the names of required inputs.

- modeldata:

  Modeldata object in the calling function.

- calling_env:

  Calling environment, should be
  [`rlang::caller_env()`](https://rlang.r-lib.org/reference/stack.html).

## Details

The difference between `tbp` (to be provided) and `tbc` (to be computed)
is that `tbp` is for promises based on data and assumptions, and `tbc`
is for promises based on internal modeldata variables (typically in
`.metainfo`). This means that there can also be cascades: once a `tbp`
is resolved, this may lead to an update of some metainformation, such
that a `tbc` can be resolved in the next update.
