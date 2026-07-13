# Computes a digest of the stan model source files, i.e. of the main stan file and all included stan files (without requiring the model to be compiled)

Computes a digest of the stan model source files, i.e. of the main stan
file and all included stan files (without requiring the model to be
compiled)

## Usage

``` r
get_stan_files_digest(
  modelfolder = "inst/stan",
  modelname = "EpiSewer_main.stan"
)
```
