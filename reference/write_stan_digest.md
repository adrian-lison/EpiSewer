# Compute a digest of the stan model and write it to a text file

Computes a digest (checksum) of the main stan model file and all
included stan files (e.g. functions), and writes it to `stan_digest.txt`
in the model folder. This digest is stored in every
[`EpiSewerJob()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewerJob.md)
created with the package, so that it can be verified against the digest
of the stan model compiled into a docker image (see
[`sewer_pull_docker()`](https://adrian-lison.github.io/EpiSewer/reference/sewer_pull_docker.md))
before fitting a job in that image.

## Usage

``` r
write_stan_digest(modelfolder = "inst/stan", modelname = "EpiSewer_main.stan")
```

## Arguments

- modelfolder:

  The folder containing the stan model file.

- modelname:

  The name of the stan model file.
