#' Summarize parameters of interest
#'
#' @description This function summarizes important parameters of interest from a
#'   fitted `EpiSewer` model.
#'
#' @param fit The fitted `EpiSewer` model.
#' @param data The model data that the `EpiSewer` model was fitted with.
#' @param .metainfo Meta information about the model data.
#' @param intervals The credible intervals (CrIs) that should be calculated. By
#'   default, these are the 50% and 95% CrIs.
#' @param ndraws Number of exemplary posterior samples that should be extracted.
#'   Note that the summaries always use all available samples.
#'
#' @details Normally, this function is automatically called by EpiSewer after
#'   model fitting, but it may also be run manually, for example to update the
#'   summary intervals or number of exemplary samples. Note that this is only
#'   possible if the result contains the full fitted model object (see
#'   [set_results_opts()]).
#'
#' @details The summaries for infections and R include a column `seeding`, which
#'   indicates whether the corresponding date was still in the seeding phase or
#'   not. The seeding phase is here defined as 2xG days long, where G is the
#'   maximum generation time. This is in contrast to the model, where it is only
#'   G days long. The reason for this decision is that after G days, new
#'   infections are still based on the seeded infections and not strongly
#'   informed by the data. By using the extended seeding criterion, we ensure
#'   that only sufficiently data-driven estimates are shown in the result plots.
#'
#' @return A `list` with the following summaries (each in a `data.frame`):
#' - R (posterior summary)
#' - R_samples (exemplary samples)
#' - expected_infections (posterior summary)
#' - infections (posterior summary)
#' - infections_samples (exemplary samples)
#' - expected_load (posterior summary)
#' - expected_concentration (posterior summary)
#' - concentration (posterior summary)
#' - sample_effects (posterior summary)
#' @export
summarize_fit <- function(fit, data, .metainfo, intervals = c(0.5, 0.95), ndraws = 50) {
  summary <- list()
  T_shift_R <- with(data, L + S + D - (G + se))
  T_shift_latent <- with(data, L + S + D)
  T_shift_onset <- with(data, S)
  T_shift_load <- with(data, 0)

  index_cols <- c(".draw", "date", "type", "seeding")

  # pseudo-randomize selected draws
  set.seed(which.min(fit$draws("R")))
  draw_ids <- sample.int(prod(dim(fit$draws("R[1]"))), ndraws)

  # draws
  R_samples <- get_latent_trajectories(
    fit, var = "R", var_forecast = "R_forecast",
    T_shift = T_shift_R,
    .metainfo = .metainfo
  )

  expected_I_samples <- get_latent_trajectories(
    fit,  var = "iota", var_forecast = "iota_forecast",
    T_shift = T_shift_latent, .metainfo = .metainfo
  )

  if (data$I_sample) {
    I_samples <- get_latent_trajectories(
      fit,  var = "I", var_forecast = "I_forecast",
      T_shift = T_shift_latent, .metainfo = .metainfo
    )
  }  else {
    I_samples <- setnames(copy(expected_I_samples), "iota", "I")
  }

  all_samples <- Reduce(
    function(...) merge(..., all = TRUE, by = c(".draw", "date", "type")),
    list(R_samples, expected_I_samples, I_samples)
    )
  min_date <- all_samples[, min(date)]
  last_seeding <- min_date + (data$G * 2) + data$se # se is seeding extension
  all_samples[, seeding := date < last_seeding]
  setcolorder(all_samples, c(".draw", "date", "type", "seeding"))
  all_samples[, type := factor(type, levels = c("estimate", "forecast"))]

  gen_dist <- data$generation_dist

  ## add growth rate
  all_samples[, infectiousness := frollapply(
      iota, FUN = weighted.mean, n = length(gen_dist),
      weights = gen_dist
      ), by = ".draw"]
  all_samples[, growth_rate := (
    log(infectiousness) - shift(log(infectiousness), type="lag")
    ), by = ".draw"]
  all_samples[, growth_rate := shift(
    growth_rate, type="lead",
    n = round(weighted.mean(1:length(gen_dist), gen_dist))
    ), by = ".draw"]

  ## add doubling time
  all_samples[ , doubling_time := log(2)/growth_rate]

  ## add days_growing
  all_samples[
    shift(R, type="lead") > 1 & R < 0, date_R_cross_1 := date - min_date
  ]
  all_samples[R <= 1, date_R_cross_1 := date - min_date]
  all_samples[
    shift(date, type="lead") == min_date + length(gen_dist) &
      shift(R, type="lead") > 1,
    date_R_cross_1 := date - min_date
    ]
  all_samples[
    , date_R_cross_1 := nafill(date_R_cross_1, type = "locf"), by = ".draw"
  ]
  all_samples[, date_R_cross_1 := min_date + date_R_cross_1]
  all_samples[, days_growing := date - date_R_cross_1]
  all_samples[, date_R_cross_1 := NULL]

  summary[["samples"]] <- all_samples[.draw %in% draw_ids, ]

  # summaries

  summary[["R"]] <- summarize_samples_dt(
    all_samples, index_cols = index_cols,
    variable = "R", intervals = intervals, cols_end = c(2,3)
  )[!(is.na(median) & type == "estimate")]

  summary[["R_diagnostics"]] <- get_diagnostics_1d_date(
    fit, var = "R", var_forecast = "R_forecast",
    T_shift = T_shift_R, .metainfo = .metainfo,
    intervals = intervals
  )
  summary[["R_diagnostics"]][, seeding := summary[["R"]][, "seeding"]]

  summary[["expected_infections"]] <- summarize_samples_dt(
    all_samples, index_cols = index_cols,
    variable = "iota", intervals = intervals, cols_end = c(2,3)
  )[!(is.na(median) & type == "estimate")]

  summary[["infections"]] <- summarize_samples_dt(
    all_samples, index_cols = index_cols,
    variable = "I", intervals = intervals, cols_end = c(2,3)
  )[!(is.na(median) & type == "estimate")]

  summary[["growth_rate"]] <- summarize_samples_dt(
    all_samples, index_cols = index_cols,
    variable = "growth_rate", intervals = intervals, cols_end = c(2,3)
  )[!(is.na(median) & type == "estimate")]

  summary[["doubling_time"]] <- summarize_samples_dt(
    all_samples, index_cols = index_cols,
    variable = "doubling_time", intervals = intervals, cols_end = c(2,3)
  )[!(is.na(median) & type == "estimate")]

  summary[["days_growing"]] <- summarize_samples_dt(
    all_samples, index_cols = index_cols,
    variable = "days_growing", intervals = intervals, cols_end = c(2,3)
  )[!(is.na(median) & type == "estimate")]

  days_growing_probs <- all_samples[
    , setNames(lapply(
      c(3,7,14,21,28),
      function(d) sum(days_growing>=d, na.rm = T)/sum(!is.na(days_growing))
      ), paste("at_least", c(3,7,14,21,28), sep = "_")),
    by = c("date", "type", "seeding")
  ]

  summary[["days_growing"]] <- merge(
    summary[["days_growing"]], days_growing_probs,
    by = c("date", "type", "seeding")
  )

  summary[["expected_load"]] <- get_summary_1d_date_log(
    fit, "pi_log",
    T_shift = T_shift_load, .metainfo = .metainfo,
    var_forecast = "pi_log_forecast",
    intervals = intervals
  )

  summary[["expected_concentration"]] <- get_summary_1d_date_log(
    fit, "kappa_log",
    T_shift = T_shift_load, .metainfo = .metainfo,
    var_forecast = "kappa_log_forecast",
    intervals = intervals
  )

  summary[["concentration"]] <- get_summary_1d_date(
    fit, "predicted_concentration",
    T_shift = T_shift_load, .metainfo = .metainfo,
    var_forecast = "predicted_concentration_forecast",
    intervals = intervals
  )

  summary[["normalized_concentration"]] <- get_summary_1d_date(
    fit, "predicted_concentration_norm",
    T_shift = T_shift_load, .metainfo = .metainfo,
    var_forecast = "predicted_concentration_forecast_norm",
    intervals = intervals
  )

  if (data$outliers) {
    summary[["outliers"]] <- get_summary_1d_date(
      fit, "epsilon",
      T_shift = T_shift_load, .metainfo = .metainfo,
      intervals = intervals
    )[, c("date", "median")]
    data.table::setnames(summary[["outliers"]], "median", "epsilon")
    summary[["outliers"]][, is_outlier := summary[["outliers"]][, epsilon] > 1]
    summary[["outliers"]] <- summary[["outliers"]][
      date %in% .metainfo$measured_dates,
      ]
    if(any(tail(summary[["outliers"]]$is_outlier,7))) {
      cli::cli_warn(paste(
        "One or several of the most recent 7 measurements could be an outlier.",
        "Real-time estimates may be influenced by these outliers.",
        "Please inspect `summary$outliers` for details and consider removing",
        "outliers manually before model fitting, or wait for more data before",
        "interpreting results."
        ))
    }
  }

  if (data$K > 0) {
    # here we exponentiate to get the multiplicative effect
    summary[["sample_effects"]] <- get_summary_vector_log(
      fit, "eta", colnames(data$X),
      intervals = intervals
    )
  }

  return(summary)
}
