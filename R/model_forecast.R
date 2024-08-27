#' Forecasting module
#'
#' @description This module function is used to specify the components of the
#'   `forecast` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'   functions (see available options below). See the documentation of the
#'   individual helper functions to adjust various settings.
#'
#' @param horizon The forecast horizon. How many days into the future should
#'   EpiSewer forecast? Note that this functionality is intended for short-term
#'   forecasts. Projections over longer horizons can be highly inaccurate.
#'   Available options: `r component_functions_("horizon")`
#'
#' @details Forecasts account for the estimated variation of transmission
#'   dynamics over time and therefore become more uncertain at longer forecast
#'   horizons. However, it is important to keep in mind that in expectation the
#'   model will project the current transmission dynamics to continue unchanged
#'   (with the exception of `R_estimate_ets`, which can project the current
#'   linear trend of Rt). This assumption can be violated by various factors
#'   such as depletion of susceptible individuals, changes in behavior, or
#'   public health interventions.
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `forecast` module.
#' @export
#' @family {module functions}
model_forecast <- function(horizon = horizon_none()) {
  verify_is_modeldata(horizon, "horizon")
  return(modeldata_combine(horizon))
}

#' Do not produce forecasts
#'
#' @description This option specifies that no forecasts should be made. Only
#'   estimates and concentration predictions until the last observed date are
#'   produced.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
horizon_none <- function(modeldata = modeldata_init()) {
  modeldata$h <- 0
  modeldata$.metainfo$forecast_horizon <- 0

  modeldata$.str$forecast[["horizon"]] <- list(
    horizon_none = c()
  )
  return(modeldata)
}

#' Specify the forecast horizon
#'
#' @description This option specifies a fixed forecast horizon in days. EpiSewer
#'   will produce concentration predictions and forecasts of all latent
#'   variables such as Rt, infections and load until the end of the forecast
#'   horizon.
#'
#' @param horizon The forecast horizon in days. If 0, no forecasts are produced.
#'   Note that this functionality is intended for short-term forecasts.
#'   Projections over longer horizons can be highly inaccurate.
#'
#' @inherit model_forecast details
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
horizon_assume <- function(horizon, modeldata = modeldata_init()) {
  modeldata$h <- horizon
  modeldata$.metainfo$forecast_horizon <- horizon

  modeldata$.str$forecast[["horizon"]] <- list(
    horizon_assume = c()
  )
  return(modeldata)
}
