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
plot_infections <- function(results, draws = FALSE, ndraws = NULL) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (draws) {
    data_to_plot <- combine_samples(results, "infections_samples", draws, ndraws)
  } else {
    data_to_plot <- combine_summaries(results, "infections")
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
plot_R <- function(results, draws = FALSE, ndraws = NULL) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (draws) {
    data_to_plot <- combine_samples(results, "R_samples", draws, ndraws)
  } else {
    data_to_plot <- combine_summaries(results, "R")
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
#' @param ww_data
#' @param include_noise
#'
#' @return
#' @export
#'
#' @examples
plot_concentration <- function(results, ww_data, include_noise = T) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  if (include_noise) {
    concentration_pred <- combine_summaries(results, "concentration")
  } else {
    concentration_pred <- combine_summaries(results, "expected_concentration")
  }

  concentration_measured <- vctrs::vec_rbind(!!!lapply(results, function(res) {
    return(ww_data[
      ww_data$date %in% res$meta_info$measured_dates,
      c("date", "concentration")
    ])
  }), .names_to = "model")
  concentration_measured$model <- forcats::fct_inorder(as.character(concentration_measured$model), ordered = T)

  plot <- ggplot(ww_data, aes(x = date)) +
    geom_line(aes(y = concentration), color = "grey", linetype = "dotted", size = 0.3) +
    geom_ribbon(data = concentration_pred, aes(ymin = lower, ymax = upper, fill = model), alpha = 0.3) +
    geom_point(aes(y = concentration), color = "grey", shape = 4) +
    geom_point(data = concentration_measured, aes(y = concentration)) +
    geom_line(data = concentration_pred, aes(y = median, color = model)) +
    theme_bw() +
    scale_x_date(expand = c(0, 0), date_breaks = "1 month", date_labels = "%b\n%Y") +
    xlab("Date") +
    ylab("Concentration [gene copies / L]") +
    coord_cartesian(xlim = as.Date(c(
      max(min(ww_data$date), min(concentration_pred$date)),
      min(max(ww_data$date), max(concentration_pred$date))
    )))
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
#' @param ww_data
#' @param include_noise
#'
#' @return
#' @export
#'
#' @examples
plot_load <- function(results, ww_data, include_noise = T) {
  if ("summary" %in% names(results)) {
    results <- list(results) # only one result object passed, wrap in list
  }
  load_pred <- combine_summaries(results, "expected_load")

  plot <- ggplot(ww_data, aes(x = date)) +
    geom_ribbon(data = load_pred, aes(ymin = lower, ymax = upper, fill = model), alpha = 0.3) +
    geom_line(data = load_pred, aes(y = median, color = model)) +
    theme_bw() +
    scale_x_date(expand = c(0, 0), date_breaks = "1 month", date_labels = "%b\n%Y") +
    xlab("Date") +
    ylab("Load [gene copies / day]") +
    coord_cartesian(xlim = as.Date(c(
      max(min(ww_data$date), min(load_pred$date)),
      min(max(ww_data$date), max(load_pred$date))
    )))
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
    geom_pointinterval(aes(x = median - 1, xmin = lower - 1, xmax = upper - 1), fatten_point = 5, ) +
    theme_bw() +
    xlab("Percentage change") +
    ylab("Predictor") +
    scale_x_continuous(labels = scales::percent) +
    scale_y_discrete(limits = rev)
}
