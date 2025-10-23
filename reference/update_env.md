# Add required data and assumptions to an environment

Add required data and assumptions to an environment

## Usage

``` r
update_env(
  f_env = list(),
  data = list(),
  assumptions = list(),
  required_data = c(),
  required_assumptions = c()
)
```

## Arguments

- f_env:

  The environment to be updated

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

The updated environment, where all required data and assumptions that
were provided in `data` and `assumptions` were inserted if missing.

## Details

Note that data and assumptions are only assigned to the environment if
the corresponding attributes are not yet present. If already present,
they are not overwritten.
