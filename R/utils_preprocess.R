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
#' To determine how much deviation from the median is ' significant, a moving
#median absolute deviation (as a more robust estimate than the standard
#deviation of how much noise to expect) in measurements is used. This ' seems to
#be more robust than just multiplying the median with a factor to ' determine
#the threshold. The moving MAD is lagged by one day such that the current value
#is ' not included. Moreover, note that because the window for the moving median
#is centered. the last window_size/2 dates have no spike detection.
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
#' @export
mark_outlier_spikes_median <- function(df, col, date_col = date, window = 5, threshold_factor = 5, mad_window = 14, mad_lower_quantile = 0.05) {
  median_info <- df %>% 
    group_by({{date_col}}) %>% summarize(daily_median = median({{col}}), .groups = "drop") %>% 
    transmute({{date_col}},
              rolling_median = zoo::rollmedian(daily_median, window, align = "center", fill = NA),
              rolling_mad = zoo::rollapply(lag(daily_median), mad_window, FUN = sd, align = "right", fill = NA),
              lower_rolling_mad = quantile(rolling_mad, mad_lower_quantile, na.rm = T)
    )
  df <- df %>% left_join(median_info, by = rlang::as_string(ensym(date_col))) %>% 
    mutate(is_outlier = {{col}} - rolling_median > threshold_factor * pmax(rolling_mad, lower_rolling_mad)) %>%
    select(-c(rolling_median, rolling_mad, lower_rolling_mad))
  return(df)
}
