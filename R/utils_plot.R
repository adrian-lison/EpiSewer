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
      aes(ymin = lower_0.95, ymax = upper_0.95, fill = model),
      alpha = 0.2, color = NA
    ) +
    geom_ribbon(
      data = data_estimate,
      aes(ymin = lower_0.5, ymax = upper_0.5, fill = model),
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
        aes(ymin = lower_0.95, ymax = upper_0.95, fill = model),
        alpha = 0.1, color = NA
      ) +
      geom_ribbon(
        data = data_forecast,
        aes(ymin = lower_0.5, ymax = upper_0.5, fill = model),
        alpha = 0.3, color = NA
      ) +
      geom_line(
        data = data_forecast,
        aes(y = upper_0.95, color = model),
        alpha = 0.4, linetype = "dashed", size = 0.2
      )  +
      geom_line(
        data = data_forecast,
        aes(y = lower_0.95, color = model),
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
      aes(ymin = lower_0.95, ymax = upper_0.95),
      alpha = 0.2, color = NA, fill = "black"
    ) +
    geom_ribbon(
      data = data_estimate,
      aes(ymin = lower_0.5, ymax = upper_0.5),
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
        aes(ymin = lower_0.95, ymax = upper_0.95),
        alpha = 0.1, color = NA, fill = "black"
      ) +
      geom_ribbon(
        data = data_forecast,
        aes(ymin = lower_0.5, ymax = upper_0.5),
        alpha = 0.4, color = NA, fill = "black"
      ) +
      geom_line(
        data = data_forecast,
        aes(y = upper_0.95),
        alpha = 0.3, color = "black", linetype = "dashed", size = 0.2
      )  +
      geom_line(
        data = data_forecast,
        aes(y = lower_0.95),
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

