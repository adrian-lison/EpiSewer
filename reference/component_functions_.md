# Get modeling functions for a component (internal function)

Get all available modeling functions for a given component. Argument
defaults are such that helpers are returned as a comma-separated list of
roxygen function links.

## Usage

``` r
component_functions_(
  component,
  collapse = "\n",
  prefix = "- [",
  suffix = "()]"
)
```

## Arguments

- component:

  A character vector with the name of the component

- collapse:

  A character used to collapse helpers into one string. Default is NULL
  (do not collapse).

- prefix:

  A prefix to add to each function name

- suffix:

  A suffix to add to each function name

## Value

A character vector with all available modeling functions for the
component. If collapse is not NULL, the helper functions are collapsed
into a single string.
