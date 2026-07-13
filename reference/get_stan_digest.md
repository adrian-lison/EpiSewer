# Get the digest of the stan model files of the installed EpiSewer package

Reads the digest written by
[`write_stan_digest()`](https://adrian-lison.github.io/EpiSewer/reference/write_stan_digest.md)
from the installed package. Returns `NULL` if no digest file is found,
e.g. because an older version of EpiSewer is installed that does not yet
provide one. This ensures backward compatibility with `EpiSewerJob`
objects created before this feature was introduced.

## Usage

``` r
get_stan_digest()
```

## Value

A single character string with the digest, or `NULL` if not found.
