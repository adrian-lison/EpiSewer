#' Model the sampling process
#'
#' @description This module function is used to specify the components of the
#'   `sampling` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'   functions (see available options below). See the documentation of the
#'   individual helper functions to adjust model priors and further settings.
#'
#' @param outliers Outliers in concentrations. `EpiSewer` can automatically
#'   identify independent spikes in the concentration time series and model them
#'   as outliers to reduce the impact on transmission dynamic estimates.
#'   Modeling options:
#'   `r component_functions_("outliers")`
#' @param sample_effects Sample (batch) effects. The pathogen concentration in a
#'   sample may be influenced by sampling-related external factors, for example
#'   the time between sampling and shipping to the lab (age-of-sample effect),
#'   or different sampling or storage methods. `EpiSewer` allows to estimate
#'   such effects using covariates that describe differences between the
#'   samples. Modeling options:
#'   `r component_functions_("sample_effects")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `sampling` module.
#' @export
#' @family {module functions}
model_sampling <- function(
    outliers = outliers_none(),
    sample_effects = sample_effects_none()
    ) {
  verify_is_modeldata(outliers, "outliers")
  verify_is_modeldata(sample_effects, "sample_effects")
  return(modeldata_combine(outliers, sample_effects))
}

#' Do not model sample effects
#'
#' @description This option does not model effects of sample covariates on the
#'   concentrations.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {sample effect models}
sample_effects_none <- function(modeldata = modeldata_init()) {
  modeldata$K <- 0
  modeldata$X <- numeric(0)
  modeldata$eta_prior <- numeric(0)
  modeldata$.init$eta <- numeric(0)

  modeldata$.str$sampling[["sample_effects"]] <- list(
    sample_effects_none = c()
  )

  return(modeldata)
}

#' Estimate weekday sample effects
#'
#' @description This option uses a log-linear regression model to estimate
#'   sample weekday effects on the concentration. Concentrations can be
#'   influenced by the time between sampling and shipping to the lab
#'   (age-of-sample effect), and if shipment follows a weekly batch scheme, the
#'   sampling weekday is a good proxy for the age at shipment.
#'
#' @param effect_prior_mu Prior (mean) on the regression coefficients.
#' @param effect_prior_sigma Prior (standard deviation) on the regression
#'   coefficients.
#'
#' @details Effects are estimated for weekdays Monday - Saturday, with Sunday as
#'   the baseline. `EpiSewer` will fit a fixed-effects log-linear model, random
#'   effects are currently not supported.
#'
#' @details The priors of this component have the following functional form:
#' - regression coefficients: `Normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {sample effect models}
sample_effects_estimate_weekday <- function(
    effect_prior_mu = 0,
    effect_prior_sigma = 1,
    modeldata = modeldata_init()) {
  modeldata <- tbc(
    "weekday_design_matrix",
    {
      weekdays <- lubridate::wday(
        seq.Date(
          modeldata$.metainfo$T_start_date,
          modeldata$.metainfo$T_end_date + modeldata$.metainfo$forecast_horizon,
          by = "1 day"
        ),
        label = TRUE
      )
      design_matrix <- model.matrix(
        ~wday,
        data.frame(wday = weekdays),
        contrasts.arg = list(wday = "contr.treatment")
      )[, -1]
      modeldata <-
        sample_effects_estimate_matrix(
          design_matrix, effect_prior_mu, effect_prior_sigma, modeldata
        )
    },
    required = c(".metainfo$T_start_date", ".metainfo$T_end_date", ".metainfo$forecast_horizon"),
    modeldata = modeldata
  )

  modeldata$.str$sampling[["sample_effects"]] <- list(
    sample_effects_estimate_weekday = c()
  )

  return(modeldata)
}

#' Estimate sample effects using a design matrix
#'
#' @description This option uses a log-linear regression model to estimate
#'   effects of sample covariates on the concentration. Concentrations can be
#'   influenced by sampling-related external factors, for example the time
#'   between sampling and shipping to the lab (age-of-sample effect), or
#'   different sampling or storage methods.
#'
#' @param effect_prior_mu Prior (mean) on the regression coefficients.
#' @param effect_prior_sigma Prior (standard deviation) on the regression
#'   coefficients.
#' @param design_matrix A design matrix with different covariates potentially
#'   influencing sample concentration. The matrix must have one row for each
#'   modeled day. See [stats::model.matrix()] for construction of design
#'   matrices.
#'
#' @details `EpiSewer` will fit a fixed-effects log-linear model, random effects
#'   are currently not supported.
#'
#' @details The priors of this component have the following functional form:
#' - regression coefficients: `Normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {sample effect models}
sample_effects_estimate_matrix <- function(
    design_matrix,
    effect_prior_mu = 0,
    effect_prior_sigma = 1,
    modeldata = modeldata_init()) {
  eta_prior <- set_prior(
    "eta", "normal",
    mu = effect_prior_mu, sigma = effect_prior_sigma
  )

  modeldata <- tbc(
    "check_design_matrix",
    {
      if (!(modeldata$T + modeldata$.metainfo$forecast_horizon == nrow(design_matrix))) {
        cli::cli_abort(
          paste(
            "Mismatch: Modeled time period has",
            modeldata$T + modeldata$.metainfo$forecast_horizon,
            "days (from earliest to latest date, including forecasts and",
            "accounting for the composite window length),",
            "but design matrix for sample date effects has",
            nrow(design_matrix),
            "rows."
          )
        )
      }
    },
    required = c("T", ".metainfo$forecast_horizon"),
    modeldata = modeldata
  )
  modeldata$K <- ncol(design_matrix)
  modeldata$X <- design_matrix

  modeldata$eta_prior <- eta_prior

  modeldata$.init$eta <- rep(0, modeldata$K)

  modeldata$.str$sampling[["sample_effects"]] <- list(
    sample_effects_estimate_matrix = c()
  )

  return(modeldata)
}

#' Do not model outliers
#'
#' @description This option does not model outliers in sampled concentrations.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {outlier models}
outliers_none <- function(modeldata = modeldata_init()) {
  modeldata$outliers <- FALSE
  modeldata$epsilon_prior <- numeric(0)
  modeldata$.init$epsilon <- numeric(0)

  modeldata$.str$sampling[["outliers"]] <- list(
    outliers_none = c()
  )

  return(modeldata)
}

#' Model outliers via an extreme value distribution
#'
#' @description This option models outliers in sampled concentrations using a
#'   generalized extreme value distribution.
#'
#' @param gev_prior_mu Prior (location) of the GEV distribution.
#' @param gev_prior_sigma Prior (scale) of the GEV distribution.
#' @param gev_prior_xi Prior (shape) of the GEV distribution.
#'
#' @details `EpiSewer` can automatically identify independent spikes in the
#'   concentration time series and model them as outliers to reduce the impact
#'   on transmission dynamic estimates. Moreover, when using this modeling
#'   option, a summary of potential outlier observations will be produced after
#'   model fitting (see `summary$outliers`).
#'
#' @details Naturally, the modeling of outliers works better retrospectively
#'   than in real-time. `EpiSewer` will produce a warning if it thinks that the
#'   most recent observations contain one or more outliers. In that case, it
#'   might make sense to manually remove the outlier observation or to wait
#'   until more data is available.
#'
#' @details Outliers are modeled as independent spikes in the load on a given
#'   sampling day. The spikes follow a generalized extreme value distribution
#'   (GEV) that is scaled by the average load per case. The default prior for
#'   the GEV is chosen such that distribution of spikes is extremely
#'   right-tailed, and the 95% quantile of spikes is below the load equivalent
#'   of 1 case, i.e. has minimal impact on estimated infection dynamics. Note
#'   that because of the properties of the GEV, the mean posterior estimates of
#'   loads and concentrations will be theoretically infinite if this modeling
#'   option is used. The median remains be well-behaved.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {outlier models}
outliers_estimate <- function(gev_prior_mu = 0.0025, gev_prior_sigma = 0.005,
                              gev_prior_xi = 2, modeldata = modeldata_init()) {
  modeldata$outliers <- TRUE
  modeldata$epsilon_prior <- set_prior(
    "epsilon", dist = "gev", mu = gev_prior_mu,
    sigma = gev_prior_sigma, xi = gev_prior_xi
    )
  modeldata$.init$epsilon <- tbe(
    rep(1e-4, modeldata$T + modeldata$D),
    required = c("T", "D")
  )

  modeldata$.str$sampling[["outliers"]] <- list(
    outliers_estimate = c()
  )

  return(modeldata)
}
