#' Specify forecasting functionality
#'
#' @param horizon The forecasting horizon in days. If 0 (default), no forecasts
#'   are produced.
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `forecast` module.
#' @export
model_forecast <- function(horizon = 0) {
  modeldata = modeldata_init()
  modeldata$h <- horizon
  modeldata$.metainfo$forecast_horizon <- horizon
  return(modeldata)
}
