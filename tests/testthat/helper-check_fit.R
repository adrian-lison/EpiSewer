check_result_valid <- function(result) {
  if (is.null(result)) {
    stop("Result is NULL.")
  }
  if (length(result) == 0) {
    stop("Result is empty.")
  }
  if ("errors" %in% names(result)) {
    stop("Errors in result.")
  }
  if (!all(c("summary", "diagnostics", "runtime") %in% names(result))) {
    stop("Result is missing summary, diagnostics or runtime.")
  }
  if (length(result$summary) == 0) {
    stop("Summary is empty.")
  }
  return(invisible(TRUE))
}

check_result_has_forecast <- function(result) {
  check_result_valid(result)
  if (!"R" %in% names(result$summary)) {
    stop("Result is missing reproduction number summary.")
  }
  if (result$job$metainfo$forecast_horizon == 0) {
    stop("Forecast horizon is zero.")
  }
  R_forecast <- result$summary$R[type == "forecast",]
  if (nrow(R_forecast) == 0) {
    stop("Result has no forecast entries.")
  }
  if (nrow(R_forecast) != result$job$metainfo$forecast_horizon) {
    stop("Result has wrong number of forecast dates")
  }
  return(invisible(TRUE))
}
