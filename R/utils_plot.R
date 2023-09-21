#' Title
#'
#' @param results
#' @param draws
#' @param ndraws
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
plot_infections <- function(results, draws = FALSE, ndraws = NULL, seeding = FALSE) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (draws) {
    data_to_plot <- combine_samples(
      results, "infections_samples", draws, ndraws
    )
  } else {
    data_to_plot <- combine_summaries(
      results, "infections"
    )
  }

  if (!seeding) {
    data_to_plot <- data_to_plot[!data_to_plot$seeding, ]
  }

  plot <- ggplot(data_to_plot, aes(x = date, color = model, fill = model)) +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month",
      date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Infections") +
    theme(legend.title = element_blank())

  if (draws) {
    plot <- plot +
      geom_line(aes(y = I, group = .draw), size = 0.1, alpha = 0.9)
  } else {
    plot <- plot +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
      geom_line(aes(y = median))
  }
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  plot
}

#' Title
#'
#' @param results
#' @param draws
#' @param ndraws
#'
#' @return
#' @export
#'
#' @examples
plot_R <- function(results, draws = FALSE, ndraws = NULL, seeding = FALSE) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (draws) {
    data_to_plot <- combine_samples(
      results, "R_samples", draws, ndraws
    )
  } else {
    data_to_plot <- combine_summaries(
      results, "R"
    )
  }

  if (!seeding) {
    data_to_plot <- data_to_plot[!data_to_plot$seeding, ]
  }

  plot <- ggplot(data_to_plot, aes(x = date, color = model, fill = model)) +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month",
      date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab(expression(R[t])) +
    theme(legend.title = element_blank()) +
    geom_hline(yintercept = 1, linetype = "dashed")

  if (draws) {
    plot <- plot +
      geom_line(aes(y = R, group = .draw), size = 0.1, alpha = 0.9)
  } else {
    plot <- plot +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
      geom_line(aes(y = median))
  }
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  plot
}

#' Title
#'
#' @param results
#' @param measurements
#' @param include_noise
#'
#' @return
#' @export
#'
#' @examples
plot_concentration <- function(results, measurements = NULL, include_noise = TRUE) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (include_noise) {
    concentration_pred <- combine_summaries(results, "concentration")
  } else {
    concentration_pred <- combine_summaries(results, "expected_concentration")
  }

  if (!is.null(measurements)) {
    concentration_measured <- vctrs::vec_rbind(!!!lapply(results, function(res) {
      return(measurements[
        measurements$date %in% res$job$meta_info$measured_dates,
        c("date", "concentration")
      ])
    }), .names_to = "model")
    concentration_measured$model <- forcats::fct_inorder(
      as.character(concentration_measured$model),
      ordered = TRUE
    )
  } else {
    concentration_measured <- NULL
  }

  if (is.null(measurements)) {
    first_date <- min(concentration_pred$date)
    last_date <- max(concentration_pred$date)
  } else {
    first_date <- max(min(measurements$date), min(concentration_pred$date))
    last_date <- min(max(measurements$date), max(concentration_pred$date))
  }

  if (!is.null(measurements)) {
    measurements <- measurements[
      !is.na(measurements$concentration),
    ]
    concentration_measured <- concentration_measured[
      !is.na(concentration_measured$concentration),
    ]
  }
  concentration_pred <- concentration_pred[
    !is.na(concentration_pred$median),
  ]

  plot <- ggplot(concentration_pred, aes(x = date)) +
    {
      if (!is.null(measurements)) {
        geom_line(
          data = measurements,
          aes(y = concentration),
          color = "grey", linetype = "dotted", size = 0.3
        )
      }
    } +
    geom_ribbon(
      data = concentration_pred,
      aes(ymin = lower, ymax = upper, fill = model), alpha = 0.3
    ) +
    {
      if (!is.null(measurements)) {
        geom_point(
          data = measurements,
          aes(y = concentration), color = "grey", shape = 4
        )
      }
    } +
    {
      if (!is.null(concentration_measured)) {
        geom_point(
          data = concentration_measured,
          aes(y = concentration)
        )
      }
    } +
    geom_line(
      data = concentration_pred,
      aes(y = median, color = model)
    ) +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month", date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Concentration [gene copies / L]") +
    coord_cartesian(xlim = as.Date(c(first_date, last_date)))

  if (length(unique(concentration_pred$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  plot
}

#' Title
#'
#' @param results
#' @param measurements
#' @param include_noise
#'
#' @return
#' @export
#'
#' @examples
plot_load <- function(results, measurements = NULL, include_noise = TRUE) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  load_pred <- combine_summaries(results, "expected_load")

  if (is.null(measurements)) {
    first_date <- min(load_pred$date)
    last_date <- max(load_pred$date)
  } else {
    first_date <- max(min(measurements$date), min(load_pred$date))
    last_date <- min(max(measurements$date), max(load_pred$date))
  }

  plot <- ggplot(load_pred, aes(x = date)) +
    geom_ribbon(
      data = load_pred,
      aes(ymin = lower, ymax = upper, fill = model), alpha = 0.3
    ) +
    geom_line(data = load_pred, aes(y = median, color = model)) +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month", date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Load [gene copies / day]") +
    coord_cartesian(xlim = as.Date(c(first_date, last_date)))

  if (length(unique(load_pred$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  plot
}

#' Title
#'
#' @param result
#'
#' @return
#' @export
#'
#' @examples
plot_sample_date_effects <- function(result) {
  ggplot(result$summary$sample_date_effects, aes(y = variable)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointinterval(
      aes(x = median - 1, xmin = lower - 1, xmax = upper - 1),
      fatten_point = 5,
    ) +
    theme_bw() +
    xlab("Percentage change") +
    ylab("Predictor") +
    scale_x_continuous(labels = scales::percent) +
    scale_y_discrete(limits = rev)
}

plot_LOD <- function(modeldata) {
  if (!all(c("LOD","LOD_sharpness") %in% names(modeldata))) {
    abort(c("The following variables must be present in model data:",
            "LOD", "LOD_sharpness"))
  }
  LOD_f <- function(x) {
    plogis((modeldata$LOD - x) / modeldata$LOD * modeldata$LOD_sharpness)
  }
  example_data <- data.frame(
    x = seq(modeldata$LOD-modeldata$LOD*0.75,
            modeldata$LOD+modeldata$LOD*0.75,length.out = 100)
    )
  example_data$y <- LOD_f(example_data$x)
  ggplot(example_data, aes(x=x, y=y)) +
    geom_vline(xintercept = modeldata$LOD, linetype = "dashed") +
    geom_line() +
    xlab("Measurement") +
    ylab("Probability of non-detection") +
    coord_cartesian(ylim = c(0,1)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw()
}
