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
#'   Available options:
#'   `r component_functions_("horizon")`
#' @param damping EpiSewer dampens the forecast of Rt so that trends in
#'   transmission will level off after some time. This prevents unrealistic
#'   extrapolation of transmission dynamics. Available options:
#'   `r component_functions_("damping")`
#'
#' @details Forecasts account for the estimated variation of transmission
#'   dynamics over time and therefore tend to become more uncertain at longer
#'   forecast horizons. However, it is important to keep in mind that depending
#'   on the Rt model used, EpiSewer will project the current transmission
#'   dynamics to continue unchanged (when using [R_estimate_rw()],
#'   [R_estimate_splines()], [R_estimate_piecewise()]) or according to a
#'   (dampened) linear trend (when using [R_estimate_ets()] and
#'   [R_estimate_changepoint_splines()]). This assumption can be violated by
#'   various factors such as depletion of susceptible individuals, changes in
#'   behavior, or public health interventions.
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `forecast` module.
#' @export
#' @family {module functions}
model_forecast <- function(horizon = horizon_none(),
                           damping = damping_assume(damping = 0.95)) {
  verify_is_modeldata(horizon, "horizon")
  verify_is_modeldata(damping, "damping")
  return(modeldata_combine(horizon, damping))
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

#' Do not dampen forecasts
#'
#' @description This option applies no damping to Rt forecasts, i.e. the trend
#'   projected by the R estimation model will be continued until the end of the
#'   forecast horizon.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
damping_none <- function(modeldata = modeldata_init()) {
  modeldata$forecast_damping <- 1
  modeldata$.str$forecast[["damping"]] <- list(
    damping_none = c()
  )
  return(modeldata)
}

#' Dampen forecasts
#'
#' @description This option dampens the forecast of Rt so that trends in
#'   transmission will level off after some time. This prevents unrealistic
#'   extrapolation of transmission dynamics.
#'
#' @param damping The forecast damping parameter. A value of 1 means no
#'   damping, a value of 0 means flat forecast. The default is 0.8, which
#'   levels off after a horizon of approximately 2 weeks.
#'
#' @details The applied damping is exponential, i.e. the trend is reduced by
#'   `damping^1` on the first forecast day, by `damping^2` on the second
#'   forecast day, and so on.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
damping_assume <- function(damping = 0.95, modeldata = modeldata_init()) {
  modeldata$forecast_damping <- damping
  modeldata$.str$forecast[["damping"]] <- list(
    damping_assume = c()
  )
  return(modeldata)
}
