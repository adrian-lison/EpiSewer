#' Plot infections
#'
#' @description Plots the estimated number of infections over time from a fitted
#'   `EpiSewer` model.
#'
#' @param results Results object returned by [EpiSewer()] after model fitting.
#' @param draws If `FALSE` (default), 50% and 95% Bayesian credible intervals
#'   are shown. If `TRUE`, exemplary posterior samples are shown in a "spaghetti
#'   plot" style.
#' @param ndraws Number of different samples to show if `draws=TRUE`.
#' @param seeding Should infections from the seeding phase be shown as well?
#'   Default is `FALSE`.
#' @param median Should the estimated median be shown, or only the credible
#'   intervals? Default is `FALSE` to avoid over-interpretation of the median.
#'
#' @return A ggplot object showing the time series of estimated infections,
#'   either with credible intervals or as "spaghetti plot". Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and add
#'   further geoms.
#' @export
#' @import ggplot2
plot_infections <- function(results, draws = FALSE, ndraws = NULL,
                            seeding = FALSE, median = FALSE) {
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

  ymin <- quantile(data_to_plot$lower_0.95, probs = 0.01)
  ymax <- quantile(data_to_plot$upper_0.95, probs = 0.99)

  plot <- ggplot(data_to_plot, aes(x = date, color = model, fill = model)) +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month",
      date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Infections") +
    theme(legend.title = element_blank()) +
    coord_cartesian(ylim = c(ymin, ymax))

  if (draws) {
    plot <- plot +
      geom_line(aes(y = I, group = .draw), size = 0.1, alpha = 0.9)
  } else {
    plot <- plot +
      geom_ribbon(
        aes(ymin = lower_0.95, ymax = upper_0.95),
        alpha = 0.2, color = NA
      ) +
      geom_ribbon(
        aes(ymin = lower_0.5, ymax = upper_0.5),
        alpha = 0.4, color = NA
      ) + {
        if (median) geom_line(aes(y = median))
      }
  }
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  plot
}

#' Plot the effective reproduction number
#'
#' @description Plots the effective reproduction number over time from a fitted
#'   `EpiSewer` model.
#'
#' @param results Results object returned by [EpiSewer()] after model fitting.
#' @param seeding Should Rt from the seeding phase be shown as well?
#'   Default is `FALSE`.
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the time series of estimated Rt,
#'   either with credible intervals or as "spaghetti plot". Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and add
#'   further geoms.
#' @export
plot_R <- function(results, draws = FALSE, ndraws = NULL,
                   seeding = FALSE, median = FALSE) {
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

  ymin <- min(0.7, quantile(data_to_plot$lower_0.95, probs = 0.01))
  ymax <- max(1.5, quantile(data_to_plot$upper_0.95, probs = 0.99))

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
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_cartesian(ylim = c(ymin, ymax))

  if (draws) {
    plot <- plot +
      geom_line(aes(y = R, group = .draw), size = 0.1, alpha = 0.9)
  } else {
    plot <- plot +
      geom_ribbon(
        aes(ymin = lower_0.95, ymax = upper_0.95),
        alpha = 0.2, color = NA
      ) +
      geom_ribbon(
        aes(ymin = lower_0.5, ymax = upper_0.5),
        alpha = 0.4, color = NA
      ) + {
        if (median) geom_line(aes(y = median))
      }
  }
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  plot
}

#' Plot predicted concentration
#'
#' @description Plots the predicted concentration over time from a fitted
#'   `EpiSewer` model.
#'
#' @param measurements A `data.frame` with observed measurements, which will be
#'   plotted alongside the predicted values. Useful to assess model fit.
#' @param include_noise If `TRUE` (default), concentrations including
#'   measurement noise are shown. If `FALSE`, only the expected concentrations
#'   are shown.
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing predicted and observed concentrations over
#'   time. Can be further manipulated using [ggplot2] functions to adjust themes
#'   and scales, and add further geoms.
#' @export
plot_concentration <- function(results, measurements = NULL,
                               include_noise = TRUE, median = FALSE) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (include_noise) {
    concentration_pred <- combine_summaries(results, "concentration")
  } else {
    concentration_pred <- combine_summaries(results, "expected_concentration")
  }

  if (!is.null(measurements)) {
    measurements_modeled <- vctrs::vec_rbind(!!!lapply(
      results, function(res) {
        return(measurements[
          measurements$date %in% res$job$metainfo$measured_dates,
          c("date", "concentration")
        ])
      }
    ), .names_to = "model")
    measurements_modeled$model <- forcats::fct_inorder(
      as.character(measurements_modeled$model),
      ordered = TRUE
    )
  } else {
    measurements_modeled <- NULL
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
      !is.na(measurements$concentration) &
        measurements$date<=last_date,
    ]
    measurements_modeled <- measurements_modeled[
      !is.na(measurements_modeled$concentration) &
        measurements_modeled$date<=last_date,
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
      aes(ymin = lower_0.95, ymax = upper_0.95, fill = model), alpha = 0.2
    ) +
    geom_ribbon(
      data = concentration_pred,
      aes(ymin = lower_0.5, ymax = upper_0.5, fill = model), alpha = 0.4
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
      if (!is.null(measurements_modeled)) {
        geom_point(
          data = measurements_modeled,
          aes(y = concentration), color = "black", shape = 4
        )
      }
    } +
    {
      if (median) {
        geom_line(
          data = concentration_pred,
          aes(y = median, color = model)
        )
      }
    } +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month", date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Concentration [gene copies / L]") +
    coord_cartesian(xlim = as.Date(c(first_date, last_date+1)))

  if (length(unique(concentration_pred$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  plot
}

#' Plot the estimated load
#'
#' @description Plots the estimated load in wastewater over time from a fitted
#'   `EpiSewer` model.
#'
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the estimated load over time. Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and add
#'   further geoms.
#'
#' @export
plot_load <- function(results, median = FALSE) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  load_pred <- combine_summaries(results, "expected_load")

  first_date <- min(load_pred$date)
  last_date <- max(load_pred$date)

  plot <- ggplot(load_pred, aes(x = date)) +
    geom_ribbon(
      data = load_pred,
      aes(ymin = lower_0.95, ymax = upper_0.95, fill = model), alpha = 0.2
    ) +
    geom_ribbon(
      data = load_pred,
      aes(ymin = lower_0.5, ymax = upper_0.5, fill = model), alpha = 0.4
    ) +
    {
      if (median) geom_line(data = load_pred, aes(y = median, color = model))
    } +
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

#' Plot estimated sample effects
#'
#' @description Plots estimated effect sizes for sample covariates with 95%
#'   credible intervals. Only works if the `EpiSewer` model
#'   included sample covariates, see e.g. [sample_effects_estimate_matrix()].
#'
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the estimated effect sizes. Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and add
#'   further geoms.
#' @export
plot_sample_effects <- function(results) {
  ggplot(results$summary$sample_effects, aes(y = variable)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::geom_pointinterval(
      aes(x = median - 1, xmin = lower_0.95 - 1, xmax = upper_0.95 - 1),
      fatten_point = 5,
    ) +
    theme_bw() +
    xlab("Percentage change") +
    ylab("Predictor") +
    scale_x_continuous(labels = scales::percent) +
    scale_y_discrete(limits = rev)
}


#' Plot limit of detection
#'
#' @description Helper function to visualize the assumed limit of detection.
#'
#' @param modeldata A `modeldata` object as returned by [LOD_assume()].
#'
#' @details This function can also be applied to `modeldata` objects which have
#'   been passed through several other modeling functions, as long as
#'   [LOD_assume()] was applied at some point.
#'
#' @return A plot showing the probability of non-detection (i.e. zero
#'   measurement) for different concentrations below and above the assumed LOD.
#'   Can be further manipulated using [ggplot2] functions to adjust themes and
#'   scales, and add further geoms.
#' @export
#'
#' @examples
#' modeldata <- LOD_assume(limit = 1e7, sharpness = 10)
#' plot_LOD(modeldata)
plot_LOD <- function(modeldata) {
  if (!all(c("LOD", "LOD_sharpness") %in% names(modeldata))) {
    rlang::abort(c(
      "The following variables must be present in model data:",
      "LOD", "LOD_sharpness"
    ))
  }
  LOD_f <- function(x) {
    plogis((modeldata$LOD - x) / modeldata$LOD * modeldata$LOD_sharpness)
  }
  example_data <- data.frame(
    x = seq(modeldata$LOD - modeldata$LOD * 0.75,
      modeldata$LOD + modeldata$LOD * 0.75,
      length.out = 100
    )
  )
  example_data$y <- LOD_f(example_data$x)
  ggplot(example_data, aes(x = x, y = y)) +
    geom_vline(xintercept = modeldata$LOD, linetype = "dashed") +
    geom_line() +
    xlab("True concentration") +
    ylab("Probability of non-detection") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_bw()
}
