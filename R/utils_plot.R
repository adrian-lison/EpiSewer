add_ribbons <- function(plot, data, median, has_forecast) {
  data_estimate <- data[type == "estimate",]
  data_forecast <- rbind(
    data[data[type == "estimate", .I[which.max(date)], by="model"]$V1],
    data[type == "forecast",]
  )
  date_forecast <- data[data[type == "estimate", .I[which.max(date)], by="model"]$V1]

  plot <- plot +
    geom_ribbon(
      data = data_estimate,
      aes(ymin = lower_outer, ymax = upper_outer, fill = model),
      alpha = 0.2, color = NA
    ) +
    geom_ribbon(
      data = data_estimate,
      aes(ymin = lower_inner, ymax = upper_inner, fill = model),
      alpha = 0.4, color = NA
    )

  if (median) {
    plot <- plot +
      geom_line(
        data = data_estimate,
        aes(y = median, color = model)
      )
  }

  if (has_forecast) {
    plot <- plot +
      geom_ribbon(
        data = data_forecast,
        aes(ymin = lower_outer, ymax = upper_outer, fill = model),
        alpha = 0.1, color = NA
      ) +
      geom_ribbon(
        data = data_forecast,
        aes(ymin = lower_inner, ymax = upper_inner, fill = model),
        alpha = 0.3, color = NA
      ) +
      geom_line(
        data = data_forecast,
        aes(y = upper_outer, color = model),
        alpha = 0.4, linetype = "dashed", size = 0.2
      )  +
      geom_line(
        data = data_forecast,
        aes(y = lower_outer, color = model),
        alpha = 0.4, linetype = "dashed", size = 0.2
      ) +
      geom_vline(
        data = date_forecast,
        aes(xintercept = date, color = model),
        linetype = "dotted"
      )

    if (median) {
      plot <- plot +
        geom_line(
          data = data_forecast,
          aes(y = median, color = model),
          alpha = 0.8, linetype = "dashed"
        )
    }
  }
  return(plot)
}

add_ribbons_base <- function(plot, data, median, has_forecast) {
  data_estimate <- data[type == "estimate",]
  data_forecast <- rbind(
              data[type == "estimate"][date == max(date, na.rm = T)],
              data[type == "forecast"]
            )
  date_forecast <- data[type == "estimate"][date == max(date)]

  plot <- plot +
    geom_ribbon(
      data = data_estimate,
      aes(ymin = lower_outer, ymax = upper_outer),
      alpha = 0.2, color = NA, fill = "black"
    ) +
    geom_ribbon(
      data = data_estimate,
      aes(ymin = lower_inner, ymax = upper_inner),
      alpha = 0.4, color = NA, fill = "black"
    )

  if (median) {
    plot <- plot +
      geom_line(
        data = data_estimate,
        aes(y = median),
        color = "black"
      )
  }

  if (has_forecast) {
    plot <- plot +
      geom_ribbon(
        data = data_forecast,
        aes(ymin = lower_outer, ymax = upper_outer),
        alpha = 0.1, color = NA, fill = "black"
      ) +
      geom_ribbon(
        data = data_forecast,
        aes(ymin = lower_inner, ymax = upper_inner),
        alpha = 0.4, color = NA, fill = "black"
      ) +
      geom_line(
        data = data_forecast,
        aes(y = upper_outer),
        alpha = 0.3, color = "black", linetype = "dashed", size = 0.2
      )  +
      geom_line(
        data = data_forecast,
        aes(y = lower_outer),
        alpha = 0.4, color = "black", linetype = "dashed", size = 0.2
      ) +
      geom_vline(
        data = date_forecast,
        aes(xintercept = date),
        color = "black", linetype = "dotted"
      )

    if (median) {
      plot <- plot +
        geom_line(
          data = data_forecast,
          aes(y = median),
          alpha = 0.8, color = "black", linetype = "dashed"
        )
    }
  }
  return(plot)
}

rename_intervals_plotting <- function(data_to_plot, intervals) {
  if (length(intervals) != 2) {
    cli::cli_abort(
      "Please specify exactly two intervals (inner and outer) for plotting."
    )
  } else if (any(intervals>1)) {
    cli::cli_abort(paste(
      "Please specify credible intervals between 0 and 1. For example,",
      "`intervals=c(0.5, 0.95)` corresponds to the 50% and 95% interval."
    ))
  } else {
    intervals = sort(intervals)
  }
  interval_columns <- c(paste0("lower_",intervals), paste0("upper_",intervals))
  if (!all(interval_columns %in% names(data_to_plot))) {
    cli::cli_abort(c(
      paste(
        "Cannot plot specified intervals:",
        paste(paste0(100*intervals,"%"), collapse = " and ")
        ),
      "*" = paste(
        "Please change the `intervals` argument for plotting, or add the",
        "desired summary intervals to the results (see `set_results_opts`)."
      )
    )
    )
  } else {
    # rename columns of intervals to inner and outer
    new_cols <- c("lower_inner","lower_outer","upper_inner","upper_outer")
    setnames(data_to_plot, old = interval_columns, new = new_cols)
  }
  return(data_to_plot)
}

#' Adaptive date breaks for result plots
#'
#' @description Uses monthly breaks for short time ranges and switches to
#'   wider break intervals (2, 3, 4, 6 months, or 1, 2, 5, 10 years) for longer
#'   time ranges, such that the number of breaks stays readable.
#'
#' @param max_breaks Maximum number of breaks to show.
#'
#' @return A function that takes the axis limits and returns date breaks.
#' @keywords internal
date_breaks_adaptive <- function(max_breaks = 8) {
  function(limits) {
    limits <- as.Date(limits)
    if (any(!is.finite(limits))) {
      return(NULL)
    }
    span_days <- as.numeric(diff(range(limits)))
    if (span_days < 60) {
      # too short for monthly breaks
      return(scales::breaks_pretty()(limits))
    }
    start <- lubridate::floor_date(min(limits), "month")
    end <- lubridate::ceiling_date(max(limits), "month")
    candidates <- c(1, 2, 3, 4, 6, 12, 24, 60, 120)
    for (k in candidates) {
      if (k < 12) {
        start_k <- lubridate::floor_date(start, paste(k, "months"))
        breaks <- seq(start_k, end, by = paste(k, "months"))
      } else {
        start_k <- lubridate::floor_date(start, paste(k / 12, "years"))
        breaks <- seq(start_k, end, by = paste(k / 12, "years"))
      }
      breaks <- breaks[breaks >= min(limits) & breaks <= max(limits)]
      if (length(breaks) <= max_breaks) {
        return(breaks)
      }
    }
    breaks
  }
}

#' Labels for adaptive date breaks
#'
#' @description Labels breaks by month and year. The year is only shown below
#'   the month for the first break and whenever the year changes. If all breaks
#'   fall on January 1st (yearly breaks), only the year is shown.
#'
#' @param breaks A vector of dates.
#'
#' @return A character vector of labels.
#' @keywords internal
date_labels_adaptive <- function(breaks) {
  breaks <- as.Date(breaks)
  labels <- rep(NA_character_, length(breaks))
  valid <- !is.na(breaks)
  if (!any(valid)) {
    return(labels)
  }
  b <- breaks[valid]
  if (all(format(b, "%d") == "01") && length(b) > 1) {
    if (all(format(b, "%m") == "01")) {
      labels[valid] <- format(b, "%Y")
    } else {
      years <- format(b, "%Y")
      new_year <- c(TRUE, years[-1] != years[-length(years)])
      labels[valid] <- ifelse(
        new_year, format(b, "%b\n%Y"), format(b, "%b")
      )
    }
  } else {
    # non-monthly breaks (short time ranges)
    labels[valid] <- format(b, "%d %b\n%Y")
  }
  labels
}

#' Adaptive date scale for result plots
#'
#' @param ... Further arguments passed to [ggplot2::scale_x_date()].
#' @keywords internal
scale_x_date_adaptive <- function(...) {
  ggplot2::scale_x_date(
    breaks = date_breaks_adaptive(),
    labels = date_labels_adaptive,
    ...
  )
}
