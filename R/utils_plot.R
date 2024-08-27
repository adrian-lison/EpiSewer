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
