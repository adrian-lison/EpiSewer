#' Compute crude estimate of load by date of infections
#'
#' @description The function smoothes the measured loads using a loess smoother
#'   on the log scale and shifts them backward by the mean total delay to obtain
#'   a crude estimate of the loads by date of infection. This is useful for
#'   different calibration purposes, such as determining a suitable
#'   `load_per_case` assumption.
#'
#' @param max_shift How far into the past should the time series be expanded? If
#'   NULL (default), only expands by the mean total delay. If `max_shift` is
#'   supplied and larger than he mean total delay, LOESS extrapolation is used.
#' @param impute_zero How should zero measurements be imputed? If NA (default),
#'   they are removed. This is the imputation value for the measured
#'   concentration. Imputed loads are calculated as `impute_zero * flow`.
#' @param impute_zero_runs Should runs of zero measurements be imputed by
#'   accounting for how many measurements in a row have been zero? This is based
#'   on the posterior expectation of the concentration for a number of zero
#'   measurements. In effect, the imputed concentration will get smaller and
#'   smaller during a run of zeros.
#' @param interpolate Should a log-linear interpolation of missing dates be
#'   performed? This makes the loess smoothing window more consistent.
#' @param loess_window Length of the loess smoothing window in days. Will be
#'   converted to the span parameter.
#' @param plot_smoothed_curve Should the smoothed curve be plotted for
#'   inspection purposes?
#'
#' @return A `data.frame` with columns `date_index`, `date`, and `load`.
#' @keywords internal
get_load_curve_crude <- function(
    measured_concentrations, measure_to_sample, sample_to_date, flow,
    total_delay_dist, max_shift = NULL, T_start_date, impute_zero = NA,
    impute_zero_runs = FALSE,
    interpolate = TRUE, loess_window = 56, plot_smoothed_curve = FALSE
    ) {
  # Estimate total delay between infections and concentration at sampling site
  total_delay <- round(dist_mean(total_delay_dist))

  # Construct a data frame with concentrations and loads
  loads <- data.table(
    date = sample_to_date[measure_to_sample],
    concentration = measured_concentrations,
    flow = flow[sample_to_date[measure_to_sample]]
  )
  loads[, load := concentration * flow]

  load_all_dates <- seq(min(loads$date, na.rm = T), max(loads$date, na.rm = T), 1)

  # Impute zero measurements
  if (is.na(impute_zero)) {
    loads <- loads[concentration > 0]
  } else if (impute_zero > 0) {
    if (impute_zero_runs) {
      n_run <- 0
      for (t in 1:nrow(loads)) {
        if (loads[t, concentration] == 0) {
          n_run <- n_run + 1
          # if n zeros have been observed in a row, the posterior expectation
          # for the concentration will be 1/(LOD_expected_scale*n)
          loads[t, concentration := impute_zero/n_run]
        } else {
            n_run <- 0
        }
      }
    } else {
      loads[load == 0, concentration := impute_zero]
    }
    # update imputed loads
    loads[, load := concentration * flow]
  } else {
    cli::cli_abort("Impute zero must be NA or a positive number.")
  }

  # Compute daily average
  loads <- loads[, .(load = mean(load, na.rm = TRUE)), by = date]

  # Log-linear interpolation of missing dates
  if (interpolate) {
    loads <- data.table(
      date = load_all_dates,
      interpolated = !load_all_dates %in% loads$dates,
      load = exp(approx(loads$date, log(loads$load), xout = load_all_dates)$y)
    )
  } else {
    loads[,interpolated := FALSE]
  }

  # Smooth loads using LOESS
  loess_fit <- loess(
    formula = log(load) ~ date,
    data = loads,
    na.action = na.exclude,
    span = loess_window / nrow(loads),
    control = loess.control(surface = "direct")
    )
  if (is.null(max_shift)) {
    max_shift <- total_delay
  }
  dates_predict <- seq(
    min(loads$date) - (max_shift - total_delay), max(loads$date) + total_delay, 1
    )
  loads_smoothed <- exp(predict(
    loess_fit, newdata = data.table(date = dates_predict), se = FALSE
    ))

  # Calculate the load curve by date of infections
  loads_by_infection <- data.table(
    date_index = dates_predict - total_delay, # accounting for delay
    date = T_start_date - 1 + dates_predict - total_delay,
    load = loads_smoothed
    )

  if (plot_smoothed_curve) {
    print(
      ggplot(loads, aes(x = date, y = load)) +
        geom_point(data = loads[interpolated==FALSE,]) +
        geom_point(data = loads[interpolated==TRUE,], color = "grey", shape = 4) +
        geom_line(
          data = data.table(date = dates_predict, load = loads_smoothed),
          color = "red"
        ) +
        geom_line(
          data = loads_by_infection, aes(x=date_index),
          color = "blue"
        ) +
        xlab("Date") + ylab("Load") +
        theme_minimal()
    )
  }

  return(loads_by_infection)
}

#' Compute crude estimate of infections by date of infections
#'
#' @description The function smoothed the measured loads using a loess smoother
#'   on the log scale and shifts them backward by the mean total delay to obtain
#'   a crude estimate of the loads by date of infection. This is useful for
#'   different calibration purposes, such as determining a suitable
#'   `load_per_case` assumption.
#'
#' @param load_curve_crude A `data.frame` with columns `date_index`, `date`, and
#'   `load` as produced by [get_load_curve_crude()].
#' @param load_per_case Load per case in gene copies.
#'
#' @return A `data.frame` with columns `date_index`, `date`, and `load`.
#' @keywords internal
get_infection_curve_crude <- function(load_curve_crude, load_per_case) {
  infection_curve_crude <- copy(load_curve_crude)
  infection_curve_crude[, infections := load / load_per_case]
  infection_curve_crude[, load := NULL]
  return(infection_curve_crude)
}
