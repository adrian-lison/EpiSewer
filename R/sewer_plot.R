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
#' @param forecast Should forecasted infections be shown? Default is true. This
#'   requires that the model was fitted with a forecast horizon, see
#'   [model_forecast()].
#' @param forecast_horizon How many days into the future should forecasts be
#'   plotted? Note that this is restricted by the forecast horizon specified
#'   during model fitting, see [horizon_assume()].
#' @param date_margin_left By how many days into the past should the plot be
#'   expanded? Can also be negative to cut off some of the earliest dates.
#' @param date_margin_right By how many days into the future should the plot be
#'   expanded? Can also be negative to cut off some of the latest dates.
#' @param facet_models Should the plot be faceted by model? Default is `FALSE`.
#' @param facet_direction How should the facetting be done? Either in different
#'   "rows" (default) or in different "columns".
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
                            median = FALSE, seeding = FALSE,
                            forecast = TRUE, forecast_horizon = NULL,
                            date_margin_left = 0, date_margin_right = 0,
                            facet_models = FALSE, facet_direction = "rows",
                            base_model = "", model_levels = NULL,
                            intervals = c(0.5, 0.95)) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (draws) {
    data_to_plot <- combine_samples(
      results, "I", draws, ndraws
    )
  } else {
    data_to_plot <- combine_summaries(
      results, "infections"
    )
    data_to_plot <- rename_intervals_plotting(data_to_plot, intervals)
  }

  if (!is.null(model_levels)) {
    data_to_plot$model <- factor(
      data_to_plot$model, levels = model_levels, ordered = TRUE
      )
  }

  if (!seeding) {
    data_to_plot <- data_to_plot[seeding==FALSE, ]
  }

  if (!forecast) {
    data_to_plot <- data_to_plot[type == "estimate",]
  } else if (!is.null(forecast_horizon)) {
    forecast_dates <- data_to_plot[
      type == "estimate",
      .(last_forecast = lubridate::as_date(max(date, na.rm = T)) + forecast_horizon),
      by = "model"]
    data_to_plot <- merge(data_to_plot, forecast_dates, by = "model")
    data_to_plot <- data_to_plot[date <= last_forecast, ]
  }

  ymin <- quantile(data_to_plot$lower_outer, probs = 0.01, na.rm = T)
  ymax <- quantile(data_to_plot$upper_outer, probs = 0.99, na.rm = T)

  xmin <- as.Date(min(data_to_plot$date, na.rm = T)) - date_margin_left
  xmax <- as.Date(max(data_to_plot$date, na.rm = T)) + date_margin_right

  has_forecast <- "forecast" %in% data_to_plot$type

  plot <- ggplot(data_to_plot[model!=base_model,],
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
  data_base_model <- data_to_plot[model==base_model,]
  data_base_model <- data_base_model[
    , setdiff(names(data_base_model), "model"), with = FALSE
    ]

  if (draws) {
    plot <- plot +
      { if (base_model!="") {
        geom_line(
          data = data_base_model[type == "estimate",],
          aes(y = I, group = .draw),
          size = 0.1, alpha = 0.9, color = "black"
        )
      }
      } +
      {
        if (base_model!="" && has_forecast) {
          geom_line(
            data = rbind(
              data_base_model[type == "estimate"][date == max(date)],
              data_base_model[type == "forecast"]
            ),
            aes(y = I, group = .draw),
            size = 0.3, alpha = 0.9, color = "black", linetype = "dashed"
          )
        }
      } +
      geom_line(
        data = data_to_plot[model!=base_model & type == "estimate",],
        aes(y = I, group = paste0(.draw, model), color = model),
        size = 0.1, alpha = 0.9
      ) +
      {
        if (has_forecast) {
          geom_line(
            data = rbind(
              data_to_plot[model!=base_model & type == "estimate",][date == max(date)],
              data_to_plot[model!=base_model & type == "forecast",]
            ),
            aes(y = I, group = paste0(.draw, model), color = model),
            size = 0.3, alpha = 0.9, linetype = "dashed"
          )
        }
      }
  } else {
    if (base_model!="") {
      plot <- add_ribbons_base(plot, data = data_base_model, median = median, has_forecast = has_forecast)
    }
    plot <- add_ribbons(plot, data = data_to_plot[model!=base_model,], median = median, has_forecast = has_forecast)
  }
  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      ggpattern::scale_pattern_fill_manual(values = "black") +
      ggpattern::scale_pattern_color_manual(values = "black") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  } else {
    plot <- plot +
      ggpattern::scale_pattern_fill_discrete() +
      ggpattern::scale_pattern_color_discrete() +
      scale_color_discrete() +
      scale_fill_discrete()
  }
  if (facet_models) {
    if (facet_direction %in% c("cols","col")) {
      plot <- plot + facet_wrap(~model, nrow = 1)
    } else if (facet_direction %in% c("rows","row")) {
      plot <- plot + facet_wrap(~model, ncol = 1)
    } else {
      cli::cli_abort(paste0(
        'Argument `facet_direction` must be either "rows" or "cols".'
      ))
    }
    plot <- plot +
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
#' @param forecast Should forecasted Rt values be shown? Default is true. This
#'   requires that the model was fitted with a forecast horizon, see
#'   [model_forecast()].
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the time series of estimated Rt, either with
#'   credible intervals or as "spaghetti plot". Can be further manipulated using
#'   [ggplot2] functions to adjust themes and scales, and to add further geoms.
#' @export
plot_R <- function(results, draws = FALSE, ndraws = NULL,
                   median = FALSE, seeding = FALSE,
                   forecast = TRUE, forecast_horizon = NULL,
                   date_margin_left = 0, date_margin_right = 0,
                   facet_models = FALSE, facet_direction = "rows",
                   base_model = "", model_levels = NULL,
                   intervals = c(0.5, 0.95)) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (draws) {
    data_to_plot <- combine_samples(
      results, "R", draws, ndraws
    )
  } else {
    data_to_plot <- combine_summaries(
      results, "R"
    )
    data_to_plot <- rename_intervals_plotting(data_to_plot, intervals)
  }

  if (!is.null(model_levels)) {
    data_to_plot$model <- factor(
      data_to_plot$model, levels = model_levels, ordered = TRUE
      )
  }

  if (!seeding) {
    data_to_plot <- data_to_plot[seeding==FALSE, ]
  }

  if (!forecast) {
    data_to_plot <- data_to_plot[type == "estimate",]
  } else if (!is.null(forecast_horizon)) {
    forecast_dates <- data_to_plot[
      type == "estimate",
      .(last_forecast = lubridate::as_date(max(date, na.rm = T)) + forecast_horizon),
      by = "model"]
    data_to_plot <- merge(data_to_plot, forecast_dates, by = "model")
    data_to_plot <- data_to_plot[date <= last_forecast, ]
  }

  ymin <- min(0.6, quantile(data_to_plot$lower_outer, probs = 0.01, na.rm = T))
  ymax <- max(1.6, quantile(data_to_plot$upper_outer, probs = 0.99, na.rm = T))

  xmin <- as.Date(min(data_to_plot$date, na.rm = T)) - date_margin_left
  xmax <- as.Date(max(data_to_plot$date, na.rm = T)) + date_margin_right

  has_forecast <- "forecast" %in% data_to_plot$type
  model_horizon <-
    rbind(
      data_to_plot[data_to_plot[type == "estimate", .I[which.max(date)], by="model"]$V1],
      data_to_plot[data_to_plot[type == "forecast", .I[which.max(date)], by="model"]$V1]
    )
  model_horizon <- model_horizon[, .(h = as.numeric(max(date) - min(date))), by = "model"]

  if (draws) {
    data_to_plot <- data_to_plot[!is.na(R),]
  } else {
    data_to_plot <- data_to_plot[!is.na(median),]
  }

  plot <- ggplot(data_to_plot[model!=base_model,],
                 aes(x = date)) +
    theme_bw() +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month",
      date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab("Effective reproduction number") +
    theme(legend.title = element_blank()) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))

  if (!base_model %in% c("", as.character(data_to_plot$model))) {
    cli::cli_abort(paste0(
      'Base model "', base_model,
      '" could not be found in the provided `results` list.'
    ))
  }
  data_base_model <- data_to_plot[model==base_model,]
  data_base_model <- data_base_model[
    , setdiff(names(data_base_model), "model"), with = FALSE
    ]

  if (draws) {
    plot <- plot +
      { if (base_model!="") {
        geom_line(
          data = data_base_model[type == "estimate",],
          aes(y = R, group = .draw),
          size = 0.1, alpha = 0.9, color = "black"
          )
        }
      } +
      {
        if (base_model!="" && has_forecast) {
          geom_line(
            data = rbind(
              data_base_model[type == "estimate"][date == max(date)],
              data_base_model[type == "forecast"]
            ),
            aes(y = R, group = .draw),
            size = 0.3, alpha = 0.9, color = "black", linetype = "dashed"
          )
        }
      } +
      geom_line(
        data = data_to_plot[model!=base_model & type == "estimate",],
        aes(y = R, group = paste0(.draw, model), color = model),
        size = 0.1, alpha = 0.9
        ) +
      {
        if (has_forecast) {
          geom_line(
            data = rbind(
              data_to_plot[model!=base_model & type == "estimate",][date == max(date)],
              data_to_plot[model!=base_model & type == "forecast",]
              ),
            aes(y = R, group = paste0(.draw, model), color = model),
            size = 0.3, alpha = 0.9, linetype = "dashed"
          )
        }
      }
  } else {
    if (base_model!="") {
      plot <- add_ribbons_base(plot, data = data_base_model, median = median, has_forecast = has_forecast)
    }
    plot <- add_ribbons(plot, data = data_to_plot[model!=base_model,], median = median, has_forecast = has_forecast)
  }

  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      ggpattern::scale_pattern_fill_manual(values = "black") +
      ggpattern::scale_pattern_color_manual(values = "black") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  } else {
    plot <- plot +
      ggpattern::scale_pattern_fill_discrete() +
      ggpattern::scale_pattern_color_discrete() +
      scale_color_discrete() +
      scale_fill_discrete()
  }
  if (facet_models) {
    if (facet_direction %in% c("cols","col")) {
      plot <- plot + facet_wrap(~model, nrow = 1)
    } else if (facet_direction %in% c("rows","row")) {
      plot <- plot + facet_wrap(~model, ncol = 1)
    } else {
      cli::cli_abort(paste0(
        'Argument `facet_direction` must be either "rows" or "cols".'
      ))
    }
    plot <- plot +
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
#' @param flows A `data.frame` with observed flows. If `normalized = TRUE`, this
#'   will be used to normalize the estimated and observed concentrations to the
#'   median flow.
#' @param include_noise If `TRUE` (default), concentrations including
#'   measurement noise are shown. If `FALSE`, only the expected concentrations
#'   are shown.
#' @param normalized If `TRUE`, the estimated and observed concentrations are
#'   normalized to the median flow. Thereby, the noise in concentrations that is
#'   due to variation in flow is essentially removed. This is especially useful
#'   for assessing forecasting performance when the future flow values were not
#'   known during model fitting.
#' @param forecast Should forecasted concentrations be shown? Default is true.
#'   This requires that the model was fitted with a forecast horizon, see
#'   [model_forecast()]. Because concentrations depend on the flow volume,
#'   accurate concentration forecasts are only possible if flow values beyond
#'   the estimation date are provided. This is typically only possible when
#'   estimating retrospectively. If flow values beyond the estimation date are
#'   not available, the forecasted concentrations will be based on the median
#'   flow volume.
#' @param mark_outliers If `TRUE`, outliers in the `measurements` are
#'   highlighted in red. See also argument `outlier_col` below.
#' @param concentration_col Name of the column in the measurements `data.frame`
#'   which contains the measured concentration that should be plotted.
#' @param date_col Name of the date column in the measurements `data.frame`.
#' @param outlier_col Name of a logical column in the measurements `data.frame`
#'   which identifies outlier measurements (for example added by
#'   [mark_outlier_spikes_median()]).
#' @param date_margin_right By how many days into the future should the plot be
#'   expanded? Can also be negative to cut off some of the latest dates. By
#'   default, this is 1 to improve the visibility of the latest date.
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
plot_concentration <- function(results = NULL, measurements = NULL, flows = NULL,
                               include_noise = TRUE, normalized = FALSE,
                               median = FALSE,
                               mark_outliers = FALSE,
                               forecast = TRUE,
                               forecast_horizon = NULL,
                               date_margin_left = 0, date_margin_right = 1,
                               facet_models = FALSE, facet_direction = "rows",
                               base_model = "", model_levels = NULL,
                               concentration_col = "concentration",
                               flow_col = "flow",
                               date_col = "date",
                               outlier_col = "is_outlier",
                               type = "time",
                               intervals = c(0.5, 0.95)
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

    if (normalized) {
      if (is.null(flows)) {
        cli::cli_abort(
          "To plot flow-normalized concentrations, please provide flow data via the `flow` argument."
        )
      } else {
        required_data_cols <- c(date_col, flow_col)
        if (!all(required_data_cols %in% names(flows))) {
          cli::cli_abort(
            c(paste(
              "The following columns must be present",
              "in the provided flow `data.frame`:",
              paste(required_data_cols, collapse = ", ")
            ),
            paste("Please adjust the `data.frame` or specify the right column",
                  "names via the `_col` arguments of this function."))
          )
        }

        flows = as.data.table(flows)[, .SD, .SDcols = required_data_cols]
        flows <- setnames(flows, old = c(date_col, flow_col), new = c("date", "flow"))

        flows[, date := lubridate::as_date(date)]

        if (any(duplicated(flows, by = "date"))) {
          flows <- unique(flows, by = c("date", "flow"))
          if (any(duplicated(flows, by = "date"))) { # if still duplicates
            cli::cli_abort(
              "Flow data is ambigious, duplicate dates with different values found."
            )
          }
        }
        median_flow <- median(flows$flow, na.rm = T)
        measurements <- merge(measurements, flows, by = "date")
        measurements[, concentration := concentration * flow / median_flow]
      }
    }
  }

  if (!is.null(results)) {
    if ("summary" %in% names(results)) {
      results <- list(results) # only one result object passed, wrap in list
    }
    if (include_noise) {
      if (normalized) {
        concentration_pred <- combine_summary_list(lapply(results, function(result){
          conc_norm <- copy(result$summary[["normalized_concentration"]])
          if (is.null(flows)) {
            median_flow <- median(result$job$data$flow, na.rm = T)
          }
          conc_norm[, mean := mean / median_flow]
          conc_norm[, median := median / median_flow]
          quantile_cols <- names(conc_norm)[
            names(conc_norm) %like% "lower_*" |
              names(conc_norm) %like% "upper_*"
          ]
          conc_norm[
            , (quantile_cols) := lapply(.SD, function(x) x / median_flow),
            .SDcols = quantile_cols
          ]
          return(conc_norm)
        }))
      } else {
        concentration_pred <- combine_summaries(results, "concentration")
      }
    } else {
      if (normalized) {
       cli::cli_abort("Normalized concentrations are only available with measurement noise.")
      }
      concentration_pred <- combine_summaries(results, "expected_concentration")
    }

    concentration_pred <- rename_intervals_plotting(concentration_pred, intervals)

    if (!is.null(model_levels)) {
      concentration_pred$model <- factor(
        concentration_pred$model, levels = model_levels, ordered = TRUE
        )
    }

    if (!is.null(measurements)) {
      measurements_modeled <- vctrs::vec_rbind(!!!lapply(
        results, function(res) {
          mes <- measurements[
            measurements$date %in% res$job$metainfo$measured_dates,
            .SD, .SDcols = data_col_names
          ]
          mes[, type := "observed"]
          if (res$job$metainfo$forecast_horizon > 0) {
            mes_forecast <- measurements[
              measurements$date %in% seq.Date(
                from = res$job$metainfo$T_end_date + 1,
                to = res$job$metainfo$T_end_date + res$job$metainfo$forecast_horizon,
                by = 1
              ),
              .SD, .SDcols = data_col_names
            ]
            mes_forecast[, type := "future"]
            mes <- rbind(mes, mes_forecast)
          }
          return(mes)
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

  if (!forecast) {
    concentration_pred <- concentration_pred[type == "estimate",]
  } else if (!is.null(forecast_horizon)) {
    forecast_dates <- concentration_pred[
      type == "estimate",
      .(last_forecast = lubridate::as_date(max(date, na.rm = T)) + forecast_horizon),
      by = "model"]
    concentration_pred <- merge(concentration_pred, forecast_dates, by = "model")
    concentration_pred <- concentration_pred[date <= last_forecast, ]
  }

  if (!is.null(concentration_pred)) {
    first_date <- as.Date(min(concentration_pred$date))
    last_date <- as.Date(max(concentration_pred$date))
  } else {
    first_date <-  as.Date(min(measurements$date))
    last_date <-  as.Date(max(measurements$date))
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
          , .SD, .SDcols=c("model","date", "obs_conc_ord")] ,
        by = c("model","date"), all.x = TRUE)
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

  has_forecast <- "forecast" %in% concentration_pred$type

  if (type == "time") {
    plot <- ggplot(data.frame(), aes(x = date))

    if (!is.null(concentration_pred)) {
      if (base_model!="") {
        plot <- add_ribbons_base(plot, data = concentration_pred_base_model, median = median, has_forecast = has_forecast)
      }
      plot <- add_ribbons(plot, data = concentration_pred[model!=base_model,], median = median, has_forecast = has_forecast)
    }

    plot <- plot +
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
            data = measurements_modeled[model==base_model & type == "observed",c("date","concentration")],
            aes(y = concentration), shape = 4, color = "black"
          )
        }
      } +
      {
        if (!is.null(measurements_modeled)) {
          geom_point(
            data = measurements_modeled[model==base_model & type == "future",c("date","concentration")],
            aes(y = concentration), shape = 8, color = "black"
          )
        }
      } +
      {
        if (!is.null(measurements_modeled)) {
          geom_point(
            data = measurements_modeled[model!=base_model & type == "observed"],
            aes(y = concentration, color = model), shape = 4
          )
        }
      } +
      {
        if (!is.null(measurements_modeled)) {
          geom_point(
            data = measurements_modeled[model!=base_model & type == "future"],
            aes(y = concentration, color = model), shape = 8
          )
        }
      } +
      theme_bw() +
      scale_x_date(
        expand = expansion(add=c(0,0.5)),
        date_breaks = "1 month", date_labels = "%b\n%Y"
      ) +
      xlab("Date") +
      ylab(ifelse(normalized,"Normalized concentration [gc / mL]","Concentration [gc / mL]")) +
      coord_cartesian(xlim = as.Date(c(first_date, last_date)))
  }

  if (type == "pp_check") {
    plot <- ggplot(data.frame(), aes(x = obs_conc_ord)) +
      {
        if (!is.null(concentration_pred) && base_model!="" && median) {
          ggdist::geom_pointinterval(
            orientation = "vertical", fatten_point = 0.5,
            interval_size_domain = c(0.5, 1),
            data = concentration_pred_base_model,
            aes(y=median, ymin = lower_outer, ymax = upper_outer), alpha = 0.2, color = "black"
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && base_model!="" && !median) {
          ggdist::geom_interval(
            orientation = "vertical",
            interval_size_domain = c(0.5, 1),
            data = concentration_pred_base_model,
            aes(ymin = lower_outer, ymax = upper_outer), alpha = 0.2, color = "black"
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && median) {
          ggdist::geom_pointinterval(
            orientation = "vertical", fatten_point = 0.5,
            interval_size_domain = c(0.5, 1),
            data = concentration_pred[concentration_pred$model!=base_model, ],
            aes(y=median, ymin = lower_outer, ymax = upper_outer, color = model), alpha = 0.2
          )
        }
      } +
      {
        if (!is.null(concentration_pred) && !median) {
          ggdist::geom_interval(
            orientation = "vertical", interval_size_domain = c(1, 26),
            data = concentration_pred[concentration_pred$model!=base_model, ],
            aes(ymin = lower_outer, ymax = upper_outer, color = model), alpha = 0.2
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
      ylab(ifelse(normalized,"Normalized concentration [gc / mL]","Concentration [gc / mL]"))
  }

  if (is.null(concentration_pred) || length(unique(concentration_pred$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      ggpattern::scale_pattern_fill_manual(values = "black") +
      ggpattern::scale_pattern_color_manual(values = "black") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  } else {
    plot <- plot +
      ggpattern::scale_pattern_fill_discrete() +
      ggpattern::scale_pattern_color_discrete() +
      scale_fill_discrete() +
      scale_color_discrete()
  }

  if (facet_models) {
    if (facet_direction %in% c("cols","col")) {
      plot <- plot + facet_wrap(~model, nrow = 1)
    } else if (facet_direction %in% c("rows","row")) {
      plot <- plot + facet_wrap(~model, ncol = 1)
    } else {
      cli::cli_abort(paste0(
        'Argument `facet_direction` must be either "rows" or "cols".'
      ))
    }
    plot <- plot +
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
#' @param forecast Should forecasted loads be shown? Default is true. This
#'   requires that the model was fitted with a forecast horizon, see
#'   [model_forecast()].
#'
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the estimated load over time. Can be further
#'   manipulated using [ggplot2] functions to adjust themes and scales, and to
#'   add further geoms.
#'
#' @export
plot_load <- function(results, median = FALSE,
                      forecast = TRUE, forecast_horizon = NULL,
                      date_margin_left = 0, date_margin_right = 0,
                      facet_models = FALSE, facet_direction = "rows",
                      base_model = "", model_levels = NULL,
                      intervals = c(0.5, 0.95)) {
  plot_time_series(
    results,
    variable = "expected_load",
    variable_name = "Load [gene copies / day]",
    yintercept = NULL,
    median = median,
    seeding = FALSE,
    forecast = forecast,
    forecast_horizon = forecast_horizon,
    date_margin_left = date_margin_left,
    date_margin_right = date_margin_right,
    facet_models = facet_models,
    facet_direction = facet_direction,
    base_model = base_model,
    model_levels = model_levels,
    intervals = intervals
  )
}

#' Plot the epidemic growth rate
#'
#' @description Plots the estimated growth rate over time from a fitted
#'   `EpiSewer` model.
#'
#' @param forecast Should forecasted growth rates be shown? Default is true.
#'   This requires that the model was fitted with a forecast horizon, see
#'   [model_forecast()].
#'
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the estimated growth rate over time. Can be
#'   further manipulated using [ggplot2] functions to adjust themes and scales,
#'   and to add further geoms.
#'
#' @export
plot_growth_rate <- function(results, median = FALSE, seeding = FALSE,
                      forecast = TRUE, forecast_horizon = NULL,
                      date_margin_left = 0, date_margin_right = 0,
                      facet_models = FALSE, facet_direction = "rows",
                      base_model = "", model_levels = NULL,
                      intervals = c(0.5, 0.95)) {
  plot_time_series(
    results,
    variable = "growth_rate",
    variable_name = "Growth rate",
    yintercept = 0,
    median = median,
    seeding = seeding,
    forecast = forecast,
    forecast_horizon = forecast_horizon,
    date_margin_left = date_margin_left,
    date_margin_right = date_margin_right,
    facet_models = facet_models,
    facet_direction = facet_direction,
    base_model = base_model,
    model_levels = model_levels,
    intervals = intervals
  )
}

#' Plot the epidemic doubling time
#'
#' @description Plots the estimated time-varying doubling time from a fitted
#'   `EpiSewer` model.
#'
#' @param forecast Should forecasted doubling times be shown? Default is true.
#'   This requires that the model was fitted with a forecast horizon, see
#'   [model_forecast()].
#'
#' @inheritParams plot_infections
#'
#' @return A ggplot object showing the estimated time-varying doubling time. Can
#'   be further manipulated using [ggplot2] functions to adjust themes and
#'   scales, and to add further geoms.
#'
#' @export
plot_doubling_time <- function(results, median = FALSE, seeding = FALSE,
                             forecast = TRUE, forecast_horizon = NULL,
                             date_margin_left = 0, date_margin_right = 0,
                             facet_models = FALSE, facet_direction = "rows",
                             base_model = "", model_levels = NULL,
                             intervals = c(0.5, 0.95)) {
  plt <- plot_time_series(
    results,
    variable = "doubling_time",
    variable_name = "Doubling time",
    yintercept = 0,
    median = median,
    seeding = seeding,
    forecast = forecast,
    forecast_horizon = forecast_horizon,
    date_margin_left = date_margin_left,
    date_margin_right = date_margin_right,
    facet_models = facet_models,
    facet_direction = facet_direction,
    base_model = base_model,
    model_levels = model_levels,
    intervals = intervals
  )
  plt <- suppressMessages(plt + coord_cartesian(ylim = c(-100, 100)))
}

plot_time_series <- function(results, variable, variable_name = variable,
                             yintercept = NULL,
                             median = FALSE, seeding = FALSE,
                      forecast = TRUE, forecast_horizon = NULL,
                      date_margin_left = 0, date_margin_right = 0,
                      facet_models = FALSE, facet_direction = "rows",
                      base_model = "", model_levels = NULL,
                      intervals = c(0.5, 0.95)) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  data_to_plot <- combine_summaries(results, variable)
  data_to_plot <- rename_intervals_plotting(data_to_plot, intervals)

  if (!seeding) {
    data_to_plot <- data_to_plot[seeding==FALSE, ]
  }

  if (!forecast) {
    data_to_plot <- data_to_plot[type == "estimate",]
  } else if (!is.null(forecast_horizon)) {
    forecast_dates <- data_to_plot[
      type == "estimate",
      .(last_forecast = lubridate::as_date(max(date, na.rm = T)) + forecast_horizon),
      by = "model"]
    data_to_plot <- merge(data_to_plot, forecast_dates, by = "model")
    data_to_plot <- data_to_plot[date <= last_forecast, ]
  }

  xmin <- min(data_to_plot$date, na.rm = T) - date_margin_left
  xmax <- max(data_to_plot$date, na.rm = T) + date_margin_right

  if (!is.null(model_levels)) {
    data_to_plot$model <- factor(
      data_to_plot$model, levels = model_levels, ordered = TRUE
    )
  }

  if (!base_model %in% c("", as.character(data_to_plot$model))) {
    cli::cli_abort(paste0(
      'Base model "', base_model,
      '" could not be found in the provided `results` list.'
    ))
  }
  data_base_model <- data_to_plot[model==base_model,]
  data_base_model <- data_base_model[
    ,setdiff(names(data_base_model), "model"), with = FALSE
  ]

  has_forecast <- "forecast" %in% data_to_plot$type

  plot <- ggplot(data_to_plot[model!=base_model,],
                 aes(x = date)) +
    {if(!is.null(yintercept)) geom_hline(yintercept = yintercept, linetype = "dashed") } +
    theme_bw() +
    theme(legend.title = element_blank()) +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "1 month", date_labels = "%b\n%Y"
    ) +
    xlab("Date") +
    ylab(variable_name) +
    coord_cartesian(xlim = c(xmin, xmax))

  if (base_model!="") {
    plot <- add_ribbons_base(plot, data = data_base_model, median = median, has_forecast = has_forecast)
  }
  plot <- add_ribbons(plot, data = data_to_plot[model!=base_model,], median = median, has_forecast = has_forecast)

  if (length(unique(data_to_plot$model)) == 1) {
    plot <- plot +
      theme(legend.position = "none") +
      ggpattern::scale_pattern_fill_manual(values = "black") +
      ggpattern::scale_pattern_color_manual(values = "black") +
      scale_color_manual(values = "black") +
      scale_fill_manual(values = "black")
  } else {
    plot <- plot +
      ggpattern::scale_pattern_fill_discrete() +
      ggpattern::scale_pattern_color_discrete() +
      scale_color_discrete() +
      scale_fill_discrete()
  }
  if (facet_models) {
    if (facet_direction %in% c("cols","col")) {
      plot <- plot + facet_wrap(~model, nrow = 1)
    } else if (facet_direction %in% c("rows","row")) {
      plot <- plot + facet_wrap(~model, ncol = 1)
    } else {
      cli::cli_abort(paste0(
        'Argument `facet_direction` must be either "rows" or "cols".'
      ))
    }
    plot <- plot +
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
plot_sample_effects <- function(results,
                                facet_models = FALSE,  facet_direction = "rows",
                                model_levels = NULL) {
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
      aes(x = median - 1, xmin = lower_outer - 1, xmax = upper_outer - 1, color = model),
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
    if (facet_direction %in% c("cols","col")) {
      plot <- plot + facet_wrap(~model, nrow = 1)
    } else if (facet_direction %in% c("rows","row")) {
      plot <- plot + facet_wrap(~model, ncol = 1)
    } else {
      cli::cli_abort(paste0(
        'Argument `facet_direction` must be either "rows" or "cols".'
      ))
    }
    plot <- plot +
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
      "{.code set_results_opts(fitted = TRUE)} when running `EpiSewer`."
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

#' Plot a growth report
#'
#' @description For a given date, the growth report shows the probability that
#'   infections have been growing for at least 3, 7, 14, 21, and 28 days,
#'   respectively. The report uses a diverging bar plot which is scaled between
#'   "very unlikely" (0% posterior probability) and "very likely" (100%
#'   posterior probability).
#'
#' @param result Results object returned by [EpiSewer()] after model fitting. In
#'   contrast to other plotting functions, this cannot be a list of multiple
#'   result objects, because growth report plots are always for a single model
#'   fit.
#' @param date The date for which the growth report should be plotted. If NULL,
#'   the most recent date for which a reliable report can be provided is
#'   automatically selected, see `partial_prob`.
#' @param partial_prob To select the most recent reliable date, we subtract a
#'   certain quantile of the shedding load distribution from the current date.
#'   For example, if `partial_prob=0.8` (default), we select the date for which
#'   80% of the shedding load of individuals infected before this date has been
#'   shed.
#'
#' @return A growth report plot showing the probability that infections have
#'   been growing for at least 3, 7, 14, 21, and 28 days, respectively. Can be
#'   further manipulated using [ggplot2] functions to adjust themes and scales,
#'   and to add further geoms.
#' @export
plot_growth_report <- function(result, date = NULL, partial_prob = 0.8) {
  if (!(class(result) == "list" && "summary" %in% names(result))) {
    cli::cli_abort(paste(
      "For prior-posterior visualization,",
      "please supply an `EpiSewer` results object."
    ))
  }
  if (is.null(date)) {
    # check that partial_prob between 0 and 1
    if (partial_prob < 0 || partial_prob > 1) {
      cli::cli_abort("The `partial_prob` argument must be between 0 and 1.")
    }
    tddist <- result$job$metainfo$total_delay_dist
    partial_delay <- which(cumsum(tddist) >= partial_prob)[1] - 1
    date_select <- result$job$metainfo$T_end_date - partial_delay
  } else {
    # check that date is in format %Y-%m-%d
    date_select <- tryCatch(lubridate::ymd(date), error = function(e) {
      cli::cli_abort("The date must be in format %Y-%m-%d.")
    })
  }
  days <- forcats::fct_inorder(paste(c(3,7,14,21,28), "days"), ordered = TRUE)
  result$summary$days_growing[date == date_select,] |>
    ggplot() +
    geom_hline(yintercept = -0.5, linetype = "solid") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0.5, linetype = "solid") +
    geom_bar(aes(x=days[1], y=at_least_3 - 0.5, fill = (at_least_3-0.5)), stat="identity") +
    geom_bar(aes(x=days[2], y=at_least_7 - 0.5, fill = (at_least_7-0.5)), stat="identity") +
    geom_bar(aes(x=days[3], y=at_least_14 - 0.5, fill = (at_least_14-0.5)), stat="identity") +
    geom_bar(aes(x=days[4], y=at_least_21 - 0.5, fill = (at_least_21-0.5)), stat="identity") +
    geom_bar(aes(x=days[5], y=at_least_28 - 0.5, fill = (at_least_28-0.5)), stat="identity") +
    scale_y_continuous(
      breaks = seq(-0.5, 0.5, by = 0.25),
      labels = c("very\nunlikely", "rather\nunlikely", "neutral", "rather\nlikely", "very\nlikely"),
      expand = expansion(add=c(0.02, 0.02))) +
    scale_fill_gradient2(low = "#006622", mid = "white", high = "#e60000", limits = c(-0.5, 0.5)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(face = "bold"),
      ) +
    coord_cartesian(ylim = c(-0.5,0.5)) +
    ggtitle(paste(
      "Probability that until", format(as.Date(date_select), "%b %d, %Y"),
      "infections have been growing for at least..."))
}
