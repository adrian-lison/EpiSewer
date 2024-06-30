#' Plot infections
#'
#' @description Plots the estimated number of infections over time from one or
#'   several fitted `EpiSewer` models.
#'
#' @param results Results object returned by [EpiSewer()] after model fitting.
#'   Can also be a named `list` with results from different model runs, in which
#'   case all results are plotted together and distinguished by colors.
#' @param draws If `FALSE` (default), 50% and 95% Bayesian credible intervals
#'   are shown. If `TRUE`, exemplary posterior samples are shown in a "spaghetti
#'   plot" style.
#' @param ndraws Number of different samples to show if `draws=TRUE`.
#' @param seeding Should infections from the seeding phase be shown as well?
#'   Default is `FALSE`.
#' @param median Should the estimated median be shown, or only the credible
#'   intervals? Default is `FALSE` to avoid over-interpretation of the median.
#' @param date_margin_left By how many days into the past should the plot be
#'   expanded? Can also be negative to cut off some of the earliest dates.
#' @param date_margin_right By how many days into the future should the plot be
#'   expanded? Can also be negative to cut off some of the latest dates.
#' @param facet_models Should the plot be faceted by model? Default is `FALSE`.
#' @param base_model Name of the base model (in the named list provided to
#'   `results`) which should be compared to the other models. This model will be
#'   plotted in black and will not be part of the legend.
#' @param model_levels A `character` vector with the names of the models to be
#'   included. The colors and legend will be ordered according to the order in
#'   `model_levels`.
#'
#' @return A ggplot object showing the time series of estimated infections,
#'   either with credible intervals or as "spaghetti plot". Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and to
#'   add further geoms.
#' @export
#' @import ggplot2
plot_infections <- function(results, draws = FALSE, ndraws = NULL,
                            seeding = FALSE, median = FALSE,
                            date_margin_left = 0, date_margin_right = 0,
                            facet_models = FALSE,
                            base_model = "", model_levels = NULL) {
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

  if (!is.null(model_levels)) {
    data_to_plot$model <- factor(data_to_plot$model, levels = model_levels, ordered = TRUE)
  }

  if (!seeding) {
    data_to_plot <- data_to_plot[!data_to_plot$seeding, ]
  }

  ymin <- quantile(data_to_plot$lower_0.95, probs = 0.01, na.rm = T)
  ymax <- quantile(data_to_plot$upper_0.95, probs = 0.99, na.rm = T)

  xmin <- as.Date(min(data_to_plot$date, na.rm = T)) - date_margin_left
  xmax <- as.Date(max(data_to_plot$date, na.rm = T)) + date_margin_right

  plot <- ggplot(data_to_plot[data_to_plot$model!=base_model,],
                 aes(x = date)) +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month",
      date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Infections") +
    theme(legend.title = element_blank()) +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))

  if (!base_model %in% c("", as.character(data_to_plot$model))) {
    cli::cli_abort(paste0(
      'Base model "', base_model,
      '" could not be found in the provided `results` list.'
    ))
  }
  data_base_model <- data_to_plot[data_to_plot$model==base_model,]
  data_base_model <- data_base_model[,setdiff(names(data_base_model), "model"), with = FALSE]

  if (draws) {
    plot <- plot +
      { if (base_model!="") {
        geom_line(data = data_base_model,
                  aes(y = I, group = .draw), size = 0.1, alpha = 0.9, color = "black")
      }
      } +
      geom_line(aes(y = I, group = paste0(.draw, model), color = model), size = 0.1, alpha = 0.9)
  } else {
    plot <- plot +
    { if (base_model!="") {
      geom_ribbon(
        data = data_base_model,
        aes(ymin = lower_0.95, ymax = upper_0.95),
        alpha = 0.2, color = NA, fill = "black"
      )
    }
    } +
      geom_ribbon(
        aes(ymin = lower_0.95, ymax = upper_0.95, fill = model),
        alpha = 0.2, color = NA
      ) +
      { if (base_model!="") {
        geom_ribbon(
          data = data_base_model,
          aes(ymin = lower_0.5, ymax = upper_0.5),
          alpha = 0.4, color = NA, fill = "black"
        )
      }
      } +
      geom_ribbon(
        aes(ymin = lower_0.5, ymax = upper_0.5, fill = model),
        alpha = 0.4, color = NA
      ) + {
        if (median && (base_model!="")) geom_line(data = data_base_model, aes(y = median), color = "black")
      } +
      {
        if (median) geom_line(aes(y = median, color = model))
      }
  }
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  if (facet_models) {
    plot <- plot +
      facet_wrap(~model, nrow = 1) +
      theme(
        strip.background = element_rect(fill = "white"),
        legend.position = "none"
        )
  }
  return(plot)
}

#' Plot the effective reproduction number
#'
#' @description Plots the effective reproduction number over time from one or
#'   several fitted `EpiSewer` models.
#'
#' @param seeding Should Rt from the seeding phase be shown as well? Default is
#'   `FALSE`.
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the time series of estimated Rt, either with
#'   credible intervals or as "spaghetti plot". Can be further manipulated using
#'   [ggplot2] functions to adjust themes and scales, and to add further geoms.
#' @export
plot_R <- function(results, draws = FALSE, ndraws = NULL,
                   seeding = FALSE, median = FALSE,
                   date_margin_left = 0, date_margin_right = 0,
                   facet_models = FALSE,
                   base_model = "", model_levels = NULL) {
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

  if (!is.null(model_levels)) {
    data_to_plot$model <- factor(data_to_plot$model, levels = model_levels, ordered = TRUE)
  }

  if (!seeding) {
    data_to_plot <- data_to_plot[!data_to_plot$seeding, ]
  }

  ymin <- min(0.6, quantile(data_to_plot$lower_0.95, probs = 0.01, na.rm = T))
  ymax <- max(1.6, quantile(data_to_plot$upper_0.95, probs = 0.99, na.rm = T))

  xmin <- as.Date(min(data_to_plot$date, na.rm = T)) - date_margin_left
  xmax <- as.Date(max(data_to_plot$date, na.rm = T)) + date_margin_right

  plot <- ggplot(data_to_plot[data_to_plot$model!=base_model,],
                 aes(x = date)) +
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
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))

  if (!base_model %in% c("", as.character(data_to_plot$model))) {
    cli::cli_abort(paste0(
      'Base model "', base_model,
      '" could not be found in the provided `results` list.'
    ))
  }
  data_base_model <- data_to_plot[data_to_plot$model==base_model,]
  data_base_model <- data_base_model[,setdiff(names(data_base_model), "model"), with = FALSE]

  if (draws) {
    plot <- plot +
      { if (base_model!="") {
        geom_line(data = data_base_model,
                  aes(y = R, group = .draw), size = 0.1, alpha = 0.9, color = "black")
        }
      } +
      geom_line(aes(y = R, group = paste0(.draw, model), color = model), size = 0.1, alpha = 0.9)
  } else {
    plot <- plot +
      { if (base_model!="") {
        geom_ribbon(
          data = data_base_model,
          aes(ymin = lower_0.95, ymax = upper_0.95),
          alpha = 0.2, color = NA, fill = "black"
        )
      }
      } +
      geom_ribbon(
        aes(ymin = lower_0.95, ymax = upper_0.95, fill = model),
        alpha = 0.2, color = NA
      ) +
      { if (base_model!="") {
        geom_ribbon(
          data = data_base_model,
          aes(ymin = lower_0.5, ymax = upper_0.5),
          alpha = 0.4, color = NA, fill = "black"
        )
      }
      } +
      geom_ribbon(
        aes(ymin = lower_0.5, ymax = upper_0.5, fill = model),
        alpha = 0.4, color = NA
      ) + {
        if (median && (base_model!="")) geom_line(data = data_base_model, aes(y = median), color = "black")
      } +
      {
        if (median) geom_line(aes(y = median, color = model))
      }
  }
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  if (facet_models) {
    plot <- plot +
      facet_wrap(~model, nrow = 1) +
      theme(
        strip.background = element_rect(fill = "white"),
        legend.position = "none"
      )
  }
  return(plot)
}

#' Plot predicted concentration
#'
#' @description Plots measured and/or predicted concentrations over time. The
#'   predicted concentrations are taken from a fitted `EpiSewer` model.
#'
#' @param measurements A `data.frame` with observed measurements, which will be
#'   plotted alongside the predicted values (useful to assess model fit). Can
#'   also be used with `results = NULL`, in which case only the observed
#'   measurements are plotted.
#' @param include_noise If `TRUE` (default), concentrations including
#'   measurement noise are shown. If `FALSE`, only the expected concentrations
#'   are shown.
#' @param mark_outliers If `TRUE`, outliers in the `measurements` are
#'   highlighted in red. See also argument `outlier_col` below.
#' @param concentration_col Name of the column in the measurements `data.frame`
#'   which contains the measured concentration that should be plotted.
#' @param date_col Name of the date column in the measurements `data.frame`.
#' @param outlier_col Name of a logical column in the measurements `data.frame`
#'   which identifies outlier measurements (for example added by
#'   [mark_outlier_spikes_median()]).
#' @param type If "time" (default), the concentration time series is plotted. If
#'   "pp_check", the observations are ordered by concentration and plotted
#'   against the predicted concentration (useful for posterior predictive
#'   checks).
#' @inheritParams plot_infections
#'
#' @details When plotting a posterior predictive check (`type="pp_check"`), each
#'   observed concentration is shown together with the 95% Credible Interval of
#'   the posterior predictive distribution for that observation. The
#'   observations are ordered by concentration, which helps to visualize bias
#'   and variance of the predictions as a function of the concentration.
#'
#' @return A ggplot object showing predicted and observed concentrations over
#'   time. Can be further manipulated using [ggplot2] functions to adjust themes
#'   and scales, and to add further geoms.
#' @export
plot_concentration <- function(results = NULL, measurements = NULL,
                               include_noise = TRUE, median = FALSE,
                               mark_outliers = FALSE,
                               date_margin_left = 0, date_margin_right = 0,
                               facet_models = FALSE,
                               base_model = "", model_levels = NULL,
                               concentration_col = "concentration",
                               date_col = "date",
                               outlier_col = "is_outlier",
                               type = "time"
                               ) {

  if (!(type %in% c("time","pp_check"))) {
    cli::cli_abort(paste0(
      'Argument `type` must be either "time" (for plotting concentrations over',
      ' time) or "pp_check" (for plotting a posterior predictive check).'
    ))
  }

  concentration_pred <- NULL
  measurements_modeled <- NULL

  if (!is.null(measurements)) {
    required_data_cols <- c(date_col, concentration_col)
    data_col_names <- c("date", "concentration")
    if (mark_outliers) {
      required_data_cols <- c(required_data_cols, outlier_col)
      data_col_names <- c(data_col_names, ".outlier")
    }

    if (!all(required_data_cols %in% names(measurements))) {
      cli::cli_abort(
        c(paste(
          "The following columns must be present",
          "in the provided measurements `data.frame`:",
          paste(required_data_cols, collapse = ", ")
        ),
        paste("Please adjust the `data.frame` or specify the right column",
              "names via the `_col` arguments of this function."))
      )
    }
    measurements = as.data.table(measurements)[, .SD, .SDcols = required_data_cols]
    setnames(measurements, old = required_data_cols, new = data_col_names)
  }

  if (!is.null(results)) {
    if ("summary" %in% names(results)) {
      results <- list(results) # only one result object passed, wrap in list
    }
    if (include_noise) {
      concentration_pred <- combine_summaries(results, "concentration")
    } else {
      concentration_pred <- combine_summaries(results, "expected_concentration")
    }

    if (!is.null(model_levels)) {
      concentration_pred$model <- factor(concentration_pred$model, levels = model_levels, ordered = TRUE)
    }

    if (!is.null(measurements)) {
      measurements_modeled <- vctrs::vec_rbind(!!!lapply(
        results, function(res) {
          return(measurements[
            measurements$date %in% res$job$metainfo$measured_dates,
            .SD, .SDcols = data_col_names
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
  }

  if (is.null(measurements)) {
    first_date <- as.Date(min(concentration_pred$date))
    last_date <- as.Date(max(concentration_pred$date))
  } else {
    if (!is.null(concentration_pred)) {
    first_date <- as.Date(max(min(measurements$date), min(concentration_pred$date)))
    last_date <-  as.Date(min(max(measurements$date), max(concentration_pred$date)))
    } else {
      first_date <-  as.Date(min(measurements$date))
      last_date <-  as.Date(max(measurements$date))
    }
  }
  first_date <- first_date - date_margin_left
  last_date <- last_date + date_margin_right

  if (!is.null(measurements)) {
    measurements <- measurements[
      !is.na(measurements$concentration) &
        measurements$date<=last_date,
    ]
  }
  if (!is.null(measurements_modeled)) {
    measurements_modeled <- measurements_modeled[
      !is.na(measurements_modeled$concentration) &
        measurements_modeled$date<=last_date,
    ]
    setorderv(measurements_modeled, cols = "concentration")
    measurements_modeled[
      , obs_conc_ord := forcats::fct_inorder(
        as.factor(concentration), ordered = TRUE
      )
    ]
    setorderv(measurements_modeled, cols = "date")
  }
  if (!is.null(concentration_pred)) {
  concentration_pred <- concentration_pred[
    !is.na(concentration_pred$median),
  ]
    if (!is.null(measurements_modeled)) {
      concentration_pred <- merge(
        concentration_pred, measurements_modeled[
          , .SD, .SDcols=c("date", "obs_conc_ord")] ,
        by = "date", all.x = TRUE)
    }
  }

  if (type == "pp_check") {
    if (is.null(measurements_modeled) || is.null(concentration_pred)) {
      cli::cli_abort(paste0(
        'If you want to plot a posterior predictive check (type="pp_check")',
        ', you must specify both the `results` and the `measurements` argument.'
      ))
    }
    concentration_pred <- concentration_pred[!is.na(obs_conc_ord),]
  }

  if (!base_model %in% c("", as.character(concentration_pred$model))) {
    cli::cli_abort(paste0(
      'Base model "', base_model,
      '" could not be found in the provided `results` list.'
      ))
  }
  concentration_pred_base_model <- concentration_pred[concentration_pred$model==base_model,]
  concentration_pred_base_model <- concentration_pred_base_model[
    ,setdiff(names(concentration_pred_base_model), "model"), with = FALSE
    ]

  if (type == "time") {
    plot <- ggplot(data.frame(), aes(x = date)) +
      {
        if (!is.null(measurements)) {
          geom_line(
            data = measurements,
            aes(y = concentration),
            color = "#a6a6a6", linetype = "dotted", linewidth = 0.3
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && base_model!="") {
          geom_ribbon(
            data = concentration_pred_base_model,
            aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.2, fill = "black"
          )
        }
      } +
      {
        if (!is.null(concentration_pred)) {
          geom_ribbon(
            data = concentration_pred[concentration_pred$model!=base_model, ],
            aes(ymin = lower_0.95, ymax = upper_0.95, fill = model), alpha = 0.2
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && base_model!="") {
          geom_ribbon(
            data = concentration_pred_base_model,
            aes(ymin = lower_0.5, ymax = upper_0.5), alpha = 0.2, fill = "black"
          )
        }
      } +
      {
        if (!is.null(concentration_pred)) {
          geom_ribbon(
            data = concentration_pred[concentration_pred$model!=base_model, ],
            aes(ymin = lower_0.5, ymax = upper_0.5, fill = model), alpha = 0.4
          )
        }
      } +
      {
        if (!is.null(measurements)) {
          geom_point(
            data = measurements,
            aes(y = concentration),
            color = ifelse(is.null(results),"black","#a6a6a6"), shape = 4
          )
        }
      } +
      {
        if (!is.null(measurements) && mark_outliers) {
          geom_point(
            data = measurements[.outlier == TRUE, ],
            aes(y = concentration),
            color = "red", shape = 4
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
        if (median && !is.null(concentration_pred) && base_model!="") {
          geom_line(
            data = concentration_pred_base_model,
            aes(y = median), color = "black"
          )
        }
      } +
      {
        if (median && !is.null(concentration_pred)) {
          geom_line(
            data = concentration_pred[concentration_pred$model!=base_model, ],
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
      ylab("Concentration [gene copies / mL]") +
      coord_cartesian(xlim = as.Date(c(first_date, last_date+1)))
  }

  if (type == "pp_check") {
    plot <- ggplot(data.frame(), aes(x = obs_conc_ord)) +
      {
        if (!is.null(concentration_pred) && base_model!="" && median) {
          ggdist::geom_pointinterval(
            orientation = "vertical", fatten_point = 0.5,
            interval_size_domain = c(0.5, 1),
            data = concentration_pred_base_model,
            aes(y=median, ymin = lower_0.95, ymax = upper_0.95), alpha = 0.2, color = "black"
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && base_model!="" && !median) {
          ggdist::geom_interval(
            orientation = "vertical",
            interval_size_domain = c(0.5, 1),
            data = concentration_pred_base_model,
            aes(ymin = lower_0.95, ymax = upper_0.95), alpha = 0.2, color = "black"
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && median) {
          ggdist::geom_pointinterval(
            orientation = "vertical", fatten_point = 0.5,
            interval_size_domain = c(0.5, 1),
            data = concentration_pred[concentration_pred$model!=base_model, ],
            aes(y=median, ymin = lower_0.95, ymax = upper_0.95, color = model), alpha = 0.2
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && !median) {
          ggdist::geom_interval(
            orientation = "vertical", interval_size_domain = c(1, 26),
            data = concentration_pred[concentration_pred$model!=base_model, ],
            aes(ymin = lower_0.95, ymax = upper_0.95, color = model), alpha = 0.2
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
        if (!is.null(measurements_modeled) && mark_outliers) {
          geom_point(
            data = measurements_modeled |> filter(.outlier),
            aes(y = concentration),
            color = "red", shape = 4
          )
        }
      } +
      theme_bw() +
      theme(
        axis.text.x = element_blank()
        ) +
      xlab("Observations (ordered by concentration)") +
      ylab("Concentration [gene copies / mL]")
  }

  if (is.null(concentration_pred) || length(unique(concentration_pred$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  if (facet_models) {
    plot <- plot +
      facet_wrap(~model, nrow = 1) +
      theme(
        strip.background = element_rect(fill = "white"),
        legend.position = "none"
      )
  }
  return(plot)
}

#' Plot the estimated load
#'
#' @description Plots the estimated load in wastewater over time from a fitted
#'   `EpiSewer` model.
#'
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the estimated load over time. Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and to
#'   add further geoms.
#'
#' @export
plot_load <- function(results, median = FALSE,
                      date_margin_left = 0, date_margin_right = 0,
                      facet_models = FALSE,
                      base_model = "", model_levels = NULL) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  load_pred <- combine_summaries(results, "expected_load")

  xmin <- min(load_pred$date, na.rm = T) - date_margin_left
  xmax <- max(load_pred$date, na.rm = T) + date_margin_right

  if (!is.null(model_levels)) {
    load_pred$model <- factor(load_pred$model, levels = model_levels, ordered = TRUE)
  }

  if (!base_model %in% c("", as.character(load_pred$model))) {
    cli::cli_abort(paste0(
      'Base model "', base_model,
      '" could not be found in the provided `results` list.'
    ))
  }
  data_base_model <- load_pred[load_pred$model==base_model,]
  data_base_model <- data_base_model[,setdiff(names(data_base_model), "model"), with = FALSE]

  plot <- ggplot(load_pred[load_pred$model!=base_model,],
                 aes(x = date)) +
    { if (base_model!="") {
      geom_ribbon(
        data = data_base_model,
        aes(ymin = lower_0.95, ymax = upper_0.95),
        alpha = 0.2, color = NA, fill = "black"
      )
    }
    } +
    geom_ribbon(
      aes(ymin = lower_0.95, ymax = upper_0.95, fill = model), alpha = 0.2
    ) +
    { if (base_model!="") {
      geom_ribbon(
        data = data_base_model,
        aes(ymin = lower_0.5, ymax = upper_0.5),
        alpha = 0.4, color = NA, fill = "black"
      )
    }
    } +
    geom_ribbon(
      aes(ymin = lower_0.5, ymax = upper_0.5, fill = model), alpha = 0.4
    ) +
    {
      if (median && (base_model!="")) geom_line(data = data_base_model, aes(y = median), color = "black")
    } +
    {
      if (median) geom_line(aes(y = median, color = model))
    } +
    theme_bw() +
    theme(legend.title = element_blank()) +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month", date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Load [gene copies / day]") +
    coord_cartesian(xlim = c(xmin, xmax))

  if (length(unique(load_pred$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  if (facet_models) {
    plot <- plot +
      facet_wrap(~model, nrow = 1) +
      theme(
        strip.background = element_rect(fill = "white"),
        legend.position = "none"
      )
  }
  return(plot)
}

#' Plot estimated sample effects
#'
#' @description Plots estimated effect sizes for sample covariates with 95%
#'   credible intervals. Only works if the `EpiSewer` model included sample
#'   covariates, see e.g. [sample_effects_estimate_matrix()].
#'
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the estimated effect sizes. Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and to
#'   add further geoms.
#' @export
plot_sample_effects <- function(results, facet_models = FALSE, model_levels = NULL) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  data_to_plot <- combine_summaries(results, "sample_effects")

  if (!is.null(model_levels)) {
    data_to_plot$model <- factor(data_to_plot$model, levels = model_levels, ordered = TRUE)
  }

  plot <- ggplot(data_to_plot[!is.na(data_to_plot$model),], aes(y = variable)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::geom_pointinterval(
      aes(x = median - 1, xmin = lower_0.95 - 1, xmax = upper_0.95 - 1, color = model),
      fatten_point = 5, position = position_dodge(width = 0.4)
    ) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    xlab("Percentage change") +
    ylab("Predictor") +
    scale_x_continuous(labels = scales::percent) +
    scale_y_discrete(limits = rev)
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  }
  if (facet_models) {
    plot <- plot +
      facet_wrap(~model, nrow = 1) +
      theme(
        strip.background = element_rect(fill = "white"),
        legend.position = "none"
      )
  }
  return(plot)
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
#'   scales, and to add further geoms.
#' @export
#'
#' @examples
#' modeldata <- LOD_assume(limit = 1e7, sharpness = 10)
#' plot_LOD(modeldata)
plot_LOD <- function(modeldata) {
  if (!all(c("LOD_model", "LOD_scale") %in% names(modeldata))) {
    cli::cli_abort(c(
      "The following variables must be present in model data:",
      "LOD_model", "LOD_scale"
    ))
  }
  LOD_f <- function(x) {
    exp(-x * modeldata$LOD_scale)
  }
  example_data <- data.frame(
    x = seq(0, -log(0.01)/modeldata$LOD_scale, length.out = 100)
  )
  example_data$y <- LOD_f(example_data$x)
  plot <- ggplot(example_data, aes(x = x, y = y)) +
    geom_vline(xintercept = -log(0.05)/modeldata$LOD_scale, linetype = "dashed") +
    geom_line() +
    xlab("True concentration") +
    ylab("Probability of non-detection") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(expand = expansion(add=c(0.01, 0.01))) +
    scale_x_continuous(expand = expansion(add=c(0, 0))) +
    theme_bw()
  return(plot)
}

#' Visually compare prior and posterior of a model parameter
#'
#' @param result Results object returned by [EpiSewer()] after model fitting. In
#'   contrast to other plotting functions, this cannot be a list of multiple
#'   result objects, because prior-posterior plots of multiple results
#'   simultaneously are currently not supported.
#' @param param_name Name of the single parameter to be plotted. This can either
#'   be the raw name of the parameter in the stan model, or it can be the
#'   semantic name from the single parameter dictionary (see details). Note that
#'   this must be a single parameter, it cannot be a parameter array. Also, only
#'   original parameters (not transformed parameters) for which a prior can be
#'   specified in `EpiSewer` are supported.
#'
#' @details The following parameters can be visualized (if in the model):
#' `r all_parameters(TRUE)`
#'
#' @return A plot showing the density of the prior (grey) and posterior (blue)
#'   for the respective parameter. Can be further manipulated using [ggplot2]
#'   functions to adjust themes and scales, and to add further geoms.
#' @export
plot_prior_posterior <- function(result, param_name) {
  if (!(class(result) == "list" && "summary" %in% names(result))) {
    cli::cli_abort(paste(
      "For prior-posterior visualization,",
      "please supply an `EpiSewer` results object."
      ))
  }
  if (!"fitted" %in% names(result)) {
    cli::cli_abort(paste(
      "Prior-posterior visualization is only possible if the fitted model is",
      "stored in the EpiSewer results object. Please use",
      "{.code set_fit_opts(fitted = TRUE)} when running `EpiSewer`."
      ))
  }
  all_params <- all_parameters()

  if (param_name %in% all_params$short_name) {
    param_i <- which(all_params$short_name == param_name)
    param_name <- all_params[param_i, "short_name"]
    raw_name <- all_params[param_i, "raw_name"]
    param_name_long <- all_params[param_i, "long_name"]
    param_scaling <- all_params[param_i, "scaling"]
    param_transf <- all_params[param_i, "transf"]
  } else if (param_name %in% all_parameters()$raw_name) {
    param_i <- which(all_params$raw_name == param_name)
    param_name <- all_params[param_i, "short_name"]
    raw_name <- all_params[param_i, "raw_name"]
    param_name_long <- all_params[param_i, "long_name"]
    param_scaling <- all_params[param_i, "scaling"]
    param_transf <- all_params[param_i, "transf"]
  } else {
    raw_name <- param_name
    param_scaling <- 1
    param_transf <- identity
  }

  prior_params <- result$job$data[[paste0(raw_name, "_prior")]]
  prior_dist_type <- stringr::str_extract(
    result$job$priors_text[[paste0(raw_name, "_prior_text")]], "(?<=dist = ).+?(?=,)"
  )

  tryCatch(
    posterior_draws <- as.vector(result$fitted$draws(raw_name)),
    error = function(e) {
      cli::cli_abort(paste0(
        "Parameter `", raw_name, "` not found in model fit."
      ), call = NULL)
    }
  )

  if (prior_dist_type == "normal") {
    prior_draws <- rnorm(2000, mean = prior_params[1], sd = prior_params[2])
  } else if (prior_dist_type == "truncated normal") {
    prior_draws <- extraDistr::rtnorm(
      2000, mean = prior_params[1], sd = prior_params[2], a = 0
      )
  } else if (prior_dist_type == "beta") {
    prior_draws <- rbeta(
      2000, shape1 = prior_params[1], shape2 = prior_params[2]
      )
  } else {
    cli::cli_abort(paste0(
      "Distribution type `", prior_dist_type, "` not supported for prior visualization."
    ))
  }

  # apply transformation and scaling
  prior_draws <- param_transf[[1]](prior_draws) * param_scaling[[1]]
  posterior_draws <- param_transf[[1]](posterior_draws) * param_scaling[[1]]

  x_lower = min(quantile(prior_draws, 0.01), quantile(posterior_draws, 0.01))
  x_upper = max(quantile(prior_draws, 0.99), quantile(posterior_draws, 0.99))

  plot <- ggplot() +
    geom_density((aes(x=prior_draws)), fill = "#a6a6a6", colour = "#808080", alpha = 0.7) +
    geom_density((aes(x=posterior_draws)), fill = "#273f76", colour = "#19294d",alpha = 0.7) +
    xlab(param_name_long) +
    ylab("Density") +
    theme_bw() +
    coord_cartesian(xlim = c(x_lower, x_upper)) +
    scale_x_continuous(expand = expansion(add = c(0,0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

  return(plot)
}
