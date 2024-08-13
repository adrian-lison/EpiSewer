#' Model the sampling process
#'
#' @description This module function is used to specify the components of the
#'  `sampling` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param sample_effects Sample (batch) effects. The pathogen concentration in
#' a sample may be influenced by sampling-related external factors, for example
#' the time between sampling and shipping to the lab (age-of-sample effect),
#' or different sampling or storage methods. `EpiSewer` allows to estimate
#' such effects using covariates that describe differences between the samples.
#' Modeling options:
#' `r component_functions_("sample_effects")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `sampling` module.
#' @export
model_sampling <- function(
    sample_effects = sample_effects_none()) {
  verify_is_modeldata(sample_effects, "sample_effects")
  return(modeldata_combine(sample_effects))
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
