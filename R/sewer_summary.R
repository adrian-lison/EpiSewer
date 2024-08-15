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
#'   (The summaries always use all draws.)
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
  T_shift_R <- with(data, L + S + D - G)
  T_shift_latent <- with(data, L + S + D)
  T_shift_onset <- with(data, S)
  T_shift_load <- with(data, 0)

  summary[["R"]] <- get_summary_1d_date(
    fit, "R",
    T_shift = T_shift_R, .metainfo = .metainfo,
    var_forecast = "R_forecast",
    intervals = intervals
  )
  summary[["R"]]$seeding <- FALSE
  summary[["R"]][1:(data$G), "seeding"] <- TRUE

  summary[["R_samples"]] <- get_latent_trajectories(
    fit, var = "R", var_forecast = "R_forecast",
    T_shift = T_shift_R,
    .metainfo = .metainfo, draw_ids = 1:ndraws
  )
  summary[["R_samples"]]$seeding <- FALSE
  summary[["R_samples"]][1:(data$G * ndraws), "seeding"] <- TRUE

  summary[["expected_infections"]] <- get_summary_1d_date(
    fit, "iota",
    T_shift = T_shift_latent, .metainfo = .metainfo,
    var_forecast = "iota_forecast",
    intervals = intervals
  )
  summary[["expected_infections"]]$seeding <- FALSE
  summary[["expected_infections"]][1:(data$G * 2), "seeding"] <- TRUE

  summary[["expected_infections_samples"]] <- get_latent_trajectories(
    fit,  var = "iota", var_forecast = "iota_forecast",
    T_shift = T_shift_latent, .metainfo = .metainfo, draw_ids = 1:ndraws
  )
  summary[["expected_infections_samples"]]$seeding <- FALSE
  summary[["expected_infections_samples"]][
    1:(data$G * 2 * ndraws), "seeding"
  ] <- TRUE

  if (data$I_sample) {
    summary[["infections"]] <- get_summary_1d_date(
      fit, "I",
      T_shift = T_shift_latent, .metainfo = .metainfo,
      var_forecast = "I_forecast",
      intervals = intervals
    )
    summary[["infections"]]$seeding <- FALSE
    summary[["infections"]][1:(data$G * 2), "seeding"] <- TRUE

    summary[["infections_samples"]] <- get_latent_trajectories(
      fit,  var = "I", var_forecast = "I_forecast",
      T_shift = T_shift_latent, .metainfo = .metainfo, draw_ids = 1:ndraws
    )
    summary[["infections_samples"]]$seeding <- FALSE
    summary[["infections_samples"]][
      1:(data$G * 2 * ndraws), "seeding"
    ] <- TRUE
  } else {
    summary[["infections"]] <- summary[["expected_infections"]]
    summary[["infections_samples"]] <- summary[["expected_infections_samples"]]
    colnames(summary[["infections_samples"]])[
      colnames(summary[["infections_samples"]]) == "iota"
      ] <- "I"
  }

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

  if (data$K > 0) {
    # here we exponentiate to get the multiplicative effect
    summary[["sample_effects"]] <- get_summary_vector_log(
      fit, "eta", colnames(data$X),
      intervals = intervals
    )
  }

  return(summary)
}
