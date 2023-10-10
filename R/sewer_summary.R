get_R_trajectories <- function(fit, T_shift, .metainfo, ndraws = 10) {
  fit_draws <- get_draws_1d_date(fit, "R", ndraws)
  date_mapping <- seq.Date(
    .metainfo$T_start_date - T_shift, .metainfo$T_end_date,
    by = "1 day"
  )
  fit_draws[, date := date_mapping[as.integer(date)]]
  fit_draws[, c(".chain", ".iteration", "variable") := NULL]
  return(fit_draws[])
}

get_I_trajectories <- function(fit, T_shift, .metainfo, ndraws = 10) {
  fit_draws <- get_draws_1d_date(fit, "I", ndraws)
  date_mapping <- seq.Date(
    .metainfo$T_start_date - T_shift, .metainfo$T_end_date,
    by = "1 day"
  )
  fit_draws[, date := date_mapping[as.integer(date)]]
  fit_draws[, c(".chain", ".iteration", "variable") := NULL]
  return(fit_draws[])
}

#' Summarize parameters of interest
#'
#' @description This function summarizes important parameters of interest from a
#'   fitted `EpiSewer` model.
#'
#' @param fit The fitted `EpiSewer` model.
#' @param data The model data that the `EpiSewer` model was fitted with.
#' @param .metainfo Meta information about the model data.
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
summarize_fit <- function(fit, data, .metainfo, ndraws = 50) {
  summary <- list()
  T_shift_R <- with(data, L + S + D - G)
  T_shift_latent <- with(data, L + S + D)
  T_shift_onset <- with(data, S)
  T_shift_load <- with(data, 0)

  summary[["R"]] <- get_summary_1d_date(
    fit, "R",
    T_shift = T_shift_R, .metainfo = .metainfo,
    intervals = c(0.5, 0.95)
  )
  summary[["R"]]$seeding <- FALSE
  summary[["R"]][1:(data$G), "seeding"] <- TRUE

  summary[["R_samples"]] <- get_R_trajectories(
    fit,
    T_shift = T_shift_R, .metainfo = .metainfo, ndraws = ndraws
  )
  summary[["R_samples"]]$seeding <- FALSE
  summary[["R_samples"]][1:(data$G * ndraws), "seeding"] <- TRUE

  summary[["expected_infections"]] <- get_summary_1d_date(
    fit, "iota",
    T_shift = T_shift_latent, .metainfo = .metainfo,
    intervals = c(0.5, 0.95)
  )
  summary[["expected_infections"]]$seeding <- FALSE
  summary[["expected_infections"]][1:(data$G * 2), "seeding"] <- TRUE

  if (data$I_sample) {
    summary[["infections"]] <- get_summary_1d_date(
      fit, "I",
      T_shift = T_shift_latent, .metainfo = .metainfo,
      intervals = c(0.5, 0.95)
    )
    summary[["infections"]]$seeding <- FALSE
    summary[["infections"]][1:(data$G * 2), "seeding"] <- TRUE

    summary[["infections_samples"]] <- get_I_trajectories(
      fit,
      T_shift = T_shift_latent, .metainfo = .metainfo, ndraws = ndraws
    )
    summary[["infections_samples"]]$seeding <- FALSE
    summary[["infections_samples"]][
      1:(data$G * 2 * ndraws), "seeding"
    ] <- TRUE
  }

  summary[["expected_load"]] <- get_summary_1d_date_log(
    fit, "kappa_log",
    T_shift = T_shift_load, .metainfo = .metainfo,
    intervals = c(0.5, 0.95)
  )

  summary[["expected_concentration"]] <- get_summary_1d_date_log(
    fit, "pi_log",
    T_shift = T_shift_load, .metainfo = .metainfo,
    intervals = c(0.5, 0.95)
  )

  summary[["concentration"]] <- get_summary_1d_date(
    fit, "predicted_concentration",
    T_shift = T_shift_load, .metainfo = .metainfo,
    intervals = c(0.5, 0.95)
  )

  if (data$K > 0) {
    # here we exponentiate to get the multiplicative effect
    summary[["sample_effects"]] <- get_summary_vector_log(
      fit, "eta", colnames(data$X),
      intervals = c(0.5, 0.95)
    )
  }

  return(summary)
}
