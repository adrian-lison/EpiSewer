#' Mark outlier spikes in a measurement time series
#'
#' Uses a simple heuristic to detect positive outliers (i.e. unusually high
#' spikes) in a time series. The approach is to compare a measurement with the
#' median of measurements in a small centered window. If the measurements before
#' and after are considerably lower than the current measurement, the median
#' will also be much lower. This is used as a criterion to determine outliers.
#'
#' @details To determine how much deviation from the median is significant, a
#'   moving median absolute deviation (as a more robust estimate than the
#'   standard deviation of how much noise to expect) in measurements is used.
#'   This seems to be more robust than just multiplying the median with a factor
#'   to determine the threshold. The moving MAD is lagged by one day such that
#'   the current value is not included. Moreover, note that because the window
#'   for the moving median is centered, the last window_size/2 dates have no
#'   spike detection.
#'
#' @details The method also allows for multiple measurements per day
#'   (replicates), where each replicate is evaluated individually. However, this
#'   currently does not give more weight to days with more replicates, i.e.
#'   ignores potential differences in measurement uncertainty.
#'
#' @param df A `data.frame` containing the time series of measurements.
#' @param measurement_col The name of the column with the measurements. Use
#'   dplyr-style env-variables, not characters.
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
mark_outlier_spikes_median <- function(
    df, measurement_col, date_col = date, window = 5,
    threshold_factor = 5, mad_window = 14, mad_lower_quantile = 0.05) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  median_info <- df %>%
    group_by({{ date_col }}) %>%
    summarize(daily_median = median({{ measurement_col }}), .groups = "drop") %>%
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
    mutate(is_outlier = {{ measurement_col }} - rolling_median > threshold_factor *
      pmax(rolling_mad, lower_rolling_mad)) %>%
    select(-c(rolling_median, rolling_mad, lower_rolling_mad))

  return(df)
}


#' Estimate load per case from wastewater data and case numbers
#'
#' @description This helper function uses a crude heuristic to infer the
#'   `load_per_case` based on the relationship between measured concentrations
#'   and case counts. The goal is to obtain a `load_per_case` assumption that is
#'   on the right order of magnitude - this will not be sufficient for accurate
#'   prevalence estimation from wastewater, but fully sufficient for monitoring
#'   trends and estimating Rt.
#'
#' @param measurements A `data.frame` with each row representing one
#'   measurement. Must have at least a column with dates and a column with
#'   concentration measurements. If multiple measurements per date are provided,
#'   their arithmetic mean is used.
#' @param cases A `data.frame` with each row representing one day. Must have at
#'   least a column with dates and a column with case numbers.
#' @param flows A `data.frame` with each row representing one day. Must have at
#'   least a column with dates and a column with flow measurements.
#' @param flow_constant Fixed flow volume, as an alternative to `flows`, if no
#'   regular flow measurements are available.
#' @param ascertainment_prop Proportion of all cases that get detected /
#'   reported. Can be used to account for underreporting of infections. Default
#'   is `ascertainment_prop=1`, meaning that 100% of infections become confirmed
#'   cases.
#' @param measurement_shift The specific timing between wastewater
#'   concentrations and case numbers depends on reporting delays and shedding
#'   profiles and is typically uncertain. This argument allows to shift the
#'   concentration and case number time series relative to each other and to
#'   average over several potential lags/leads, as specified by an integer
#'   vector. The default is `measurement_shift = seq(-7,7)`, i.e. a shift of
#'   concentrations between up to one week before and after case numbers.
#' @param shift_weights Weights for the shifted comparisons. Must be an numeric
#'   vector of the same length as `measurement_shift`. If `NULL` (default), the
#'   weights are chosen to be approximately inversely proportional to the shift
#'   distance.
#' @param date_col Name of the date column in all provided data frames.
#' @param concentration_col Name of the column containing the measured
#'   concentrations.
#' @param case_col Name of the column containing the case numbers.
#' @param flow_col Name of the column containing the flows.
#' @param signif_fig Significant figures to round to. Since this heuristic only
#'   provides crude estimates which should not be overinterpreted, the result
#'   gets rounded. Default is rounding to the 2 most significant figures.
#'
#' @details In the `EpiSewer` model, the `load_per_case` serves as a scaling
#'   factor describing how many pathogen particles are shed by the average
#'   infected individual overall and how much of this is detectable at the
#'   sampling site. This depends both on biological factors as well as on the
#'   specific sewage system. It is therefore almost always necessary to assume
#'   the load per case based on a comparison of measured concentrations/loads
#'   and case numbers.
#'
#' @details The heuristic used here is to fit a linear regression model with
#'   loads (computed using concentrations and flows) as dependent variable and
#'   case numbers as independent variable over all measurements. This does not
#'   explicitly account for shedding profiles or reporting delays, but the
#'   `measurement_shift` argument allows to average over a set of relative
#'   shifts between the two time series.
#'
#' @details The flow volume unit should be the same as for the concentration
#'   measurements, e.g. if concentrations are measured in gc/mL, then the flow
#'   should be in mL as well.
#'
#' @return A suggested `load_per_case` that can be used as an assumption in
#'   [load_per_case_assume()].
#' @export
#' @import data.table
suggest_load_per_case <- function(measurements, cases,
                                  flows = NULL, flow_constant = NULL,
                                  ascertainment_prop = 1,
                                  measurement_shift = seq(-7,7),
                                  shift_weights = 1/(abs(measurement_shift)+1),
                                  date_col = "date",
                                  concentration_col = "concentration",
                                  flow_col = "flow",
                                  case_col = "cases",
                                  signif_fig = 2) {

  required_data_cols <- c(date_col, concentration_col)
  if (!all(required_data_cols %in% names(measurements))) {
    rlang::abort(
      paste(
        "The following columns must be present",
        "in the provided measurements `data.frame`:",
        paste(required_data_cols, collapse = ", ")
      )
    )
  }
  measurements = as.data.table(measurements)[, .SD, .SDcols = required_data_cols]
  setnames(measurements, old = required_data_cols, new = c("date", "concentration"))
  measurements[, date := as.Date(date)]
  measurements <- measurements[
    , .(concentration = mean(concentration, na.rm=T)), by = date
    ]

  required_data_cols <- c(date_col, case_col)
  if (!all(required_data_cols %in% names(cases))) {
    rlang::abort(
      paste(
        "The following columns must be present",
        "in the provided cases `data.frame`:",
        paste(required_data_cols, collapse = ", ")
      )
    )
  }
  cases = as.data.table(cases)[, .SD, .SDcols = required_data_cols]
  setnames(cases, old = required_data_cols, new = c("date", "case_numbers"))
  cases[, date := as.Date(date)]
  cases <- cases[
    , .(case_numbers = mean(case_numbers, na.rm=T)), by = date
  ]

  if (!(ascertainment_prop>0 && ascertainment_prop<=1)) {
    rlang::abort(paste(
      "The ascertainment proportion must be between 0 and 1,",
      "e.g. 0.9 for 90%."
    ))
  }

  if (is.null(flows) && is.null(flow_constant)) {
    rlang::abort(c(
      "Please supply one of the following arguments:",
      "flows: a data frame with flow measurements",
      "flow_constant: a constant flow value"
      ))
  } else if (is.null(flows)) {
    measurements[, flow := flow_constant]
  } else {
    if (!is.null(flow_constant)) {
      rlang::warn(paste(
        "You provided both a data frame with flow measurements (flows)",
        "and a constant flow value (flow_constant).",
        "Only the `flows` argument will be used."
        ))
    }
    required_data_cols <- c(date_col, flow_col)
    if (!all(required_data_cols %in% names(flows))) {
      rlang::abort(
        paste(
          "The following columns must be present",
          "in the provided flows `data.frame`:",
          paste(required_data_cols, collapse = ", ")
        )
      )
    }
    flows = as.data.table(flows)[, .SD, .SDcols = required_data_cols]
    setnames(flows, old = required_data_cols, new = c("date", "flow"))
    flows[, date := as.Date(date)]
    flows <- flows[
      , .(flow = mean(flow, na.rm=T)), by = date
    ]
    measurements <- merge(measurements, flows, by = "date")
  }

  measurements[, load := concentration*flow]
  cases[, case_numbers := case_numbers/ascertainment_prop]
  shift_weights <- shift_weights/sum(shift_weights)

  measurements <- rbindlist(mapply(function(s, w){
    measurements[, .(
      date = date,
      load = shift(load, s, NA, "lead"),
      shift = s,
      weight = w)]
  }, s = measurement_shift, w = shift_weights, SIMPLIFY = FALSE))

  measurements <- measurements[!is.na(date) & !is.na(load)]
  cases <- cases[!is.na(date) & !is.na(case_numbers)]

  combined <- merge(measurements, cases, by = "date")
  lm_res <- lm(load ~ 0 + case_numbers, combined, weights = combined$weight)
  load_per_case <- unname(lm_res$coefficients["case_numbers"])
  # round to a certain number of significant figures
  load_per_case <- signif(load_per_case, signif_fig)
  return(load_per_case)
}
