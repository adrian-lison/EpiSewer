# Copy an environment

Copy an environment

## Usage

``` r
copy_env(env, exclude = c())
```

## Arguments

- env:

  The environment to be copied

- exclude:

  Character vector with names to be excluded in the copy

## Value

A copy of the environment

## Details

The new environment has the calling environment of `copy_env` as its
parent. This ensure that all functions in the namespace are accessible.
