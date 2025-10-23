# Check if all required data and assumptions have been provided.

Check if all required data and assumptions have been provided.

## Usage

``` r
check_provided(
  f_env = list(),
  data = list(),
  assumptions = list(),
  required_data = c(),
  required_assumptions = c()
)
```

## Arguments

- f_env:

  The environment to be checked

- data:

  A list with data

- assumptions:

  A list with assumptions

- required_data:

  A character vector with the names of required data. Different options
  for one item can be provided using "\|" as a separator. EpiSewer will
  then check if at least one of the options is present.

- required_assumptions:

  A character vector with the names of required assumptions. Different
  options for one item can be provided using "\|" as a separator.
  EpiSewer will then check if at least one of the options is present.

## Value

TRUE if all required data and assumptions are present in `f_env`.
