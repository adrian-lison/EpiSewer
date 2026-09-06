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
#' @return A module object containing the components of the
#'   `forecast` module.
#' @export
#' @family module functions
model_forecast <- function(horizon = horizon_none(),
                           damping = damping_assume(damping = 0.95)) {
  verify_is_component(horizon, "horizon")
  verify_is_component(damping, "damping")
  return(new_module("forecast", list(horizon = horizon, damping = damping)))
}

#' Do not produce forecasts
#'
#' @description This option specifies that no forecasts should be made. Only
#'   estimates and concentration predictions until the last observed date are
#'   produced.
#'
#' @inherit template_model_helpers return
#' @export
horizon_none <- function() {
  model_component("horizon_none", "horizon", {
    modeldata$h <- 0
    modeldata$.metainfo$forecast_horizon <- 0
    return(modeldata)
  })
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
#' @inherit template_model_helpers return
#' @export
horizon_assume <- function(horizon) {
  model_component("horizon_assume", "horizon", {
    modeldata$h <- horizon
    modeldata$.metainfo$forecast_horizon <- horizon
    return(modeldata)
  })
}

#' Do not dampen forecasts
#'
#' @description This option applies no damping to Rt forecasts, i.e. the trend
#'   projected by the R estimation model will be continued until the end of the
#'   forecast horizon.
#'
#' @inherit template_model_helpers return
#' @export
damping_none <- function() {
  model_component("damping_none", "damping", {
    modeldata$forecast_damping <- 1
    return(modeldata)
  })
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
#' @inherit template_model_helpers return
#' @export
damping_assume <- function(damping = 0.95) {
  model_component("damping_assume", "damping", {
    modeldata$forecast_damping <- damping
    return(modeldata)
  })
}
