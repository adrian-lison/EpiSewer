# Compile EpiSewer models

The stan models used by EpiSewer need to be compiled for your device.
This is only necessary once, after installing or updating the package.
This function compiles all models from the package.

## Usage

``` r
sewer_compile(model = NULL, force_recompile = FALSE, verbose = FALSE)
```

## Arguments

- model:

  A character vector with a specific model to compile. If NULL
  (default), all models are compiled.

- force_recompile:

  If TRUE, then models will be recompiled even if they have already been
  successfully compiled.

- verbose:

  If TRUE, warnings and detailed errors from the compilation are
  printed. This can help to diagnose compilation issues.

## Details

If one or several models are not successfully compiled, please ensure
that `cmdstan` is properly set up and try updating it to a newer version
using
[`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html).
If the problem persists, please run `sewer_compile(verbose = TRUE)` and
post the output in a new issue on GitHub, along with your
[`cmdstanr::cmdstan_version()`](https://mc-stan.org/cmdstanr/reference/set_cmdstan_path.html).
