# Check modeldata object for required variables

`modeldata_check` accepts various additional arguments which are turned
into a sophisticated error message about missing variables

## Usage

``` r
modeldata_check(
  modeldata,
  required,
  required_values = NULL,
  descriptions = modeldata_descriptions(),
  run_before = modeldata_var_requirements(),
  advice = NULL,
  throw_error = TRUE,
  calling_env = rlang::caller_env()
)
```

## Arguments

- modeldata:

  A `list` with all data variables (including priors) to be passed on to
  stan, alongside additional meta information and descriptions.

- required:

  Character vector of required variables

- descriptions:

  Optional, a list with names corresponding to all or some of the
  required variables, and values corresponding to additional
  descriptions of the variables

- run_before:

  Optional, either a character vector, or a `list` with names
  corresponding to all of the required variables, and values
  corresponding to functions that will add these variables to modeldata.

- advice:

  Additional message with advice on how to solve error.

- throw_error:

  Should an error be thrown (with verbatim error message), or should
  `TRUE`/`FALSE` be returned. Default is `TRUE`.

## Value

If throw_error is `TRUE`, nothing is returned, otherwise a logical is
returned (`TRUE` if all required variables are present)
