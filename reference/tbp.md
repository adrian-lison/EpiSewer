# Lazy-execute a function once all data and assumptions are provided

Lazy-execute a function once all data and assumptions are provided

## Usage

``` r
tbp(
  f_name,
  f_expr,
  required_data = c(),
  required_assumptions = c(),
  modeldata = NULL,
  calling_env = rlang::caller_env()
)
```

## Arguments

- f_name:

  Name of the "functions" that will be executed once everything is
  provided.

- f_expr:

  Expression with arbitrary R code in which attributes are assigned to
  `modeldata.`

- required_data:

  A character vector with the names of required data. Different options
  for one item can be provided using "\|" as a separator. EpiSewer will
  then check if at least one of the options is present.

- required_assumptions:

  A character vector with the names of required assumptions. Different
  options for one item can be provided using "\|" as a separator.
  EpiSewer will then check if at least one of the options is present.

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
