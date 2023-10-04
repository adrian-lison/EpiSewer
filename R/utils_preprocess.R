#' Mark outlier spikes in a measurement time series
#'
#' Uses a simple heuristic to detect positive outliers (i.e. unusually high
#' spikes) in a time series. The approach is to compare a measurement with the
#' median of measurements in a small centered window. If the measurements before
#' and after are considerably lower than the current measurement, the median
#' will also be much lower. This is used as a criterion to determine outliers.
#' Note that the method also allows for multiple measurements per day
#' (replicates), where each replicate is evaluated individually. However, this
#' currently does not give more weight to days with more replicates, i.e.
#' ignores potential differences in measurement uncertainty.
#'
#' To determine how much deviation from the median is significant, a moving
# median absolute deviation (as a more robust estimate than the standard
# deviation of how much noise to expect) in measurements is used. This ' seems to
# be more robust than just multiplying the median with a factor to ' determine
# the threshold. The moving MAD is lagged by one day such that the current value
# is ' not included. Moreover, note that because the window for the moving median
# is centered. the last window_size/2 dates have no spike detection.
#'
#' @param df A `data.frame` containing the time series of measurements.
#' @param col The name of the column with the measurements. Use dplyr-style
#'   env-variables, not characters.
#' @param date_col The name of the column with corresponding dates. Use
#'   dplyr-style env-variables, not characters.
#' @param window The size of the centered window for computing the median.
#' @param threshold_factor Beyond how many median absolute deviations from the
#'   median should a measurement be marked as outlier?
#' @param mad_window The size of the right-aligned, one-day-lagged window for
#'   the mean absolute deviation. This should be longer than the window for the
#'   median.
#' @param mad_lower_quantile At what quantile should the lower bound for the
#'   expected noise be? This is used to avoid false positives when concentration
#'   levels are very low.
#'
#' @return The provided `data.frame`, with an additional logical column
#'   `is_outlier`.
#'
mark_outlier_spikes_median <- function(
    df, col, date_col = date, window = 5,
    threshold_factor = 5, mad_window = 14, mad_lower_quantile = 0.05) {

  median_info <- df %>%
    group_by({{ date_col }}) %>%
    summarize(daily_median = median({{ col }}), .groups = "drop") %>%
    transmute({{ date_col }},
      rolling_median = zoo::rollmedian(
        daily_median, window,
        align = "center", fill = NA
      ),
      rolling_mad = zoo::rollapply(
        lag(daily_median), mad_window,
        FUN = sd, align = "right", fill = NA
      ),
      lower_rolling_mad = quantile(
        rolling_mad, mad_lower_quantile,
        na.rm = TRUE
      )
    )

  df <- df %>%
    left_join(median_info, by = rlang::as_string(ensym(date_col))) %>%
    mutate(is_outlier = {{ col }} - rolling_median > threshold_factor *
      pmax(rolling_mad, lower_rolling_mad)) %>%
    select(-c(rolling_median, rolling_mad, lower_rolling_mad))

  return(df)
}


#' Title
#'
#' @param measurements
#' @param flows
#' @param cases
#' @param ascertainment_prop
#' @param measurement_shift
#' @param shift_weights
#' @param date_col
#' @param measurement_col
#' @param case_col
#' @param flow_col
#' @param flow_constant
#'
#' @return
#' @export
#' @import data.table
#'
#' @examples
suggest_load_per_case <- function(measurements, cases, flows = NULL,
                                  flow_constant = NULL,
                                  ascertainment_prop = 1,
                                  measurement_shift = seq(-7,7),
                                  shift_weights = 1/(abs(measurement_shift)+1),
                                  date_col = "date",
                                  measurement_col = "concentration",
                                  flow_col = "flow",
                                  case_col = "cases") {

  required_data_cols <- c(date_col, measurement_col)
  if (!all(required_data_cols %in% names(measurements))) {
    abort(
      paste(
        "The following columns must be present",
        "in the provided measurements `data.frame`:",
        paste(required_data_cols, collapse = ", ")
      )
    )
  }
  measurements = as.data.table(measurements)
  measurements[, (date_col) := as.Date(get(date_col))]

  required_data_cols <- c(date_col, case_col)
  if (!all(required_data_cols %in% names(cases))) {
    abort(
      paste(
        "The following columns must be present",
        "in the provided cases `data.frame`:",
        paste(required_data_cols, collapse = ", ")
      )
    )
  }
  cases = as.data.table(cases)
  cases[, (date_col) := as.Date(get(date_col))]

  if (!(ascertainment_prop>0 && ascertainment_prop<=1)) {
    abort(paste(
      "The ascertainment proportion must be between 0 and 1,",
      "e.g. 0.9 for 90%."
    ))
  }

  if (is.null(flows) && is.null(flow_constant)) {
    abort(c(
      "Please supply one of the following arguments:",
      "flows: a data frame with flow measurements",
      "flow_constant: a constant flow value"
      ))
  } else if (is.null(flows)) {
    measurements[[flow_col]] <- flow_constant
  } else {
    if (!is.null(flow_constant)) {
      warn(paste(
        "You provided both a data frame with flow measurements (flows)",
        "and a constant flow value (flow_constant).",
        "Only the `flows` argument will be used."
        ))
    }
    required_data_cols <- c(date_col, flow_col)
    if (!all(required_data_cols %in% names(flows))) {
      abort(
        paste(
          "The following columns must be present",
          "in the provided measurements `data.frame`:",
          paste(required_data_cols, collapse = ", ")
        )
      )
    }
    flows = as.data.table(flows)
    flows[, (date_col) := as.Date(get(date_col))]
    measurements <- merge(measurements, flows, by = "date")
  }

  measurements[, load := get(measurement_col)*get(flow_col)]
  cases[, (case_col) := get(case_col)/ascertainment_prop]
  shift_weights <- shift_weights/sum(shift_weights)

  measurements <- rbindlist(mapply(function(s, w){
    measurements[, .(
      date = date,
      load = shift(load, s, NA, "lead"),
      shift = s,
      weight = w)]
  }, s = measurement_shift, w = shift_weights, SIMPLIFY = FALSE))

  measurements <- measurements[!is.na(get(date_col)) & !is.na(load)]
  cases <- cases[!is.na(get(date_col)) & !is.na(get(case_col))]

  combined <- merge(measurements, cases, by = "date")
  lm_res <- lm(load ~ 0 + cases, combined, weights = combined$weight)
  load_per_case <- unname(lm_res$coefficients["cases"])
  load_per_case <- signif(load_per_case, 2) # round to 2 significant figures
  return(load_per_case)
}
