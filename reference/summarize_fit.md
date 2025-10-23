# Summarize parameters of interest

This function summarizes important parameters of interest from a fitted
`EpiSewer` model.

## Usage

``` r
summarize_fit(fit, data, .metainfo, intervals = c(0.5, 0.95), ndraws = 50)
```

## Arguments

- fit:

  The fitted `EpiSewer` model.

- data:

  The model data that the `EpiSewer` model was fitted with.

- .metainfo:

  Meta information about the model data.

- intervals:

  The credible intervals (CrIs) that should be calculated. By default,
  these are the 50% and 95% CrIs.

- ndraws:

  Number of exemplary posterior samples that should be extracted. Note
  that the summaries always use all available samples.

## Value

A `list` with the following summaries (each in a `data.frame`):

- R (posterior summary)

- R_samples (exemplary samples)

- expected_infections (posterior summary)

- infections (posterior summary)

- infections_samples (exemplary samples)

- expected_load (posterior summary)

- expected_concentration (posterior summary)

- concentration (posterior summary)

- sample_effects (posterior summary)

## Details

Normally, this function is automatically called by EpiSewer after model
fitting, but it may also be run manually, for example to update the
summary intervals or number of exemplary samples. Note that this is only
possible if the result contains the full fitted model object (see
[`set_results_opts()`](https://adrian-lison.github.io/EpiSewer/reference/set_results_opts.md)).

The summaries for infections and R include a column `seeding`, which
indicates whether the corresponding date was still in the seeding phase
or not. The seeding phase is here defined as 2xG days long, where G is
the maximum generation time. This is in contrast to the model, where it
is only G days long. The reason for this decision is that after G days,
new infections are still based on the seeded infections and not strongly
informed by the data. By using the extended seeding criterion, we ensure
that only sufficiently data-driven estimates are shown in the result
plots.
