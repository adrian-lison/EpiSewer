# Plot estimated sample effects

Plots estimated effect sizes for sample covariates with 95% credible
intervals. Only works if the `EpiSewer` model included sample
covariates, see e.g.
[`sample_effects_estimate_matrix()`](https://adrian-lison.github.io/EpiSewer/reference/sample_effects_estimate_matrix.md).

## Usage

``` r
plot_sample_effects(
  results,
  facet_models = FALSE,
  facet_direction = "rows",
  model_levels = NULL
)
```

## Arguments

- results:

  Results object returned by
  [`EpiSewer()`](https://adrian-lison.github.io/EpiSewer/reference/EpiSewer.md)
  after model fitting. Can also be a named `list` with results from
  different model runs, in which case all results are plotted together
  and distinguished by colors.

- facet_models:

  Should the plot be faceted by model? Default is `FALSE`.

- facet_direction:

  How should the facetting be done? Either in different "rows" (default)
  or in different "columns".

- model_levels:

  A `character` vector with the names of the models to be included. The
  colors and legend will be ordered according to the order in
  `model_levels`.

## Value

A `ggplot` object showing the estimated effect sizes. Can be further
manipulated using `ggplot2` functions to adjust themes and scales, and
to add further geoms.
