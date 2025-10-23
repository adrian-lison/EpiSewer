# Assume a sewer residence time distribution

This option assumes a fixed residence time distribution for pathogen
particles. By default, `EpiSewer` assumes that particles arrive at the
sampling site within the day of shedding. However, for larger sewage
systems, particles may travel longer than a day depending on where and
when they were shed into the wastewater.

## Usage

``` r
residence_dist_assume(residence_dist = NULL, modeldata = modeldata_init())
```

## Arguments

- residence_dist:

  A numeric vector representing a discrete residence time distribution,
  with elements describing the share of load that takes 0 days, 1 day, 2
  days, and so on to arrive at the sampling site.

- modeldata:

  A `modeldata` object to which the above model specifications should be
  added. Default is an empty model given by
  [`modeldata_init()`](https://adrian-lison.github.io/EpiSewer/reference/modeldata_init.md).
  Can also be an already partly specified model returned by other
  `EpiSewer` modeling functions.

## Value

A `modeldata` object containing data and specifications of the model to
be fitted. Can be passed on to other `EpiSewer` modeling functions to
add further data and model specifications.

The `modeldata` object also includes information about parameter
initialization (`.init`), meta data (`.metainfo`), and checks to be
performed before model fitting (`.checks`).

## Examples

``` r
# Particles arrive within the same day
residence_dist_assume(residence_dist = c(1))
#> sewage
#>  |- residence_dist_assume

# Particles always arrive after one day
residence_dist_assume(residence_dist = c(0, 1))
#> sewage
#>  |- residence_dist_assume

# 1/4 of particles only arrives after one day
residence_dist_assume(residence_dist = c(0.75, 0.25))
#> sewage
#>  |- residence_dist_assume
```
