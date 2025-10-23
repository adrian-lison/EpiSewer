# Model outliers via an extreme value distribution

This option models outliers in sampled concentrations using a
generalized extreme value distribution.

## Usage

``` r
outliers_estimate(
  gev_prior_mu = 0,
  gev_prior_sigma = 2e-08,
  gev_prior_xi = 4,
  modeldata = modeldata_init()
)
```

## Arguments

- gev_prior_mu:

  Prior (location) of the GEV distribution.

- gev_prior_sigma:

  Prior (scale) of the GEV distribution.

- gev_prior_xi:

  Prior (shape) of the GEV distribution.

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

## Details

`EpiSewer` can automatically identify independent spikes in the
concentration time series and model them as outliers to reduce the
impact on transmission dynamic estimates. Moreover, when using this
modeling option, a summary of potential outlier observations will be
produced after model fitting (see `summary$outliers`).

Naturally, the modeling of outliers works better retrospectively than in
real-time. `EpiSewer` will produce a warning if it thinks that the most
recent observations contain one or more outliers. In that case, it might
make sense to manually remove the outlier observation or to wait until
more data is available.

Outliers are modeled as independent spikes in the load on a given
sampling day. The spikes follow a generalized extreme value distribution
(GEV) that is scaled by the average load per case. The default prior for
the GEV is chosen such that distribution of spikes is extremely
right-tailed, and the 99% quantile of spikes is below the load
equivalent of 1 case, i.e. has minimal impact on estimated infection
dynamics. Note that because of the properties of the GEV, the mean
posterior estimates of loads and concentrations will be theoretically
infinite if this modeling option is used. The median remains be
well-behaved.

## See also

Other outlier models:
[`outliers_none()`](https://adrian-lison.github.io/EpiSewer/reference/outliers_none.md)
