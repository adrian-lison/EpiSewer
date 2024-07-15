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
#' EpiSewer forecast?
#' Available options:
#' `r component_functions_("horizon")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `forecast` module.
#' @export
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
  return(horizon_assume(horizon = 0))
}

#' Specify the forecast horizon
#'
#' @description This option specifies a fixed forecast horizon in days. EpiSewer
#'   will produce concentration predictions and forecasts of all latent
#'   variables such as Rt, infections and load until the end of the forecast
#'   horizon.
#'
#' @param horizon The forecast horizon in days. If 0, no forecasts are produced.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
horizon_assume <- function(horizon, modeldata = modeldata_init()) {
  modeldata$h <- horizon
  modeldata$.metainfo$forecast_horizon <- horizon
  return(modeldata)
}
