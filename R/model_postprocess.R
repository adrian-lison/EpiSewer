get_R_trajectories <- function(fit, T_shift, meta_info, ndraws = 10) {
  fit_draws <- get_draws_1d_date(fit, "R", ndraws)
  date_mapping <- seq.Date(
    meta_info$T_start_date - T_shift, meta_info$T_end_date,
    by = "1 day"
  )
  fit_draws[, date := date_mapping[as.integer(date)]]
  fit_draws[, c(".chain", ".iteration", "variable") := NULL]
  return(fit_draws[])
}

get_I_trajectories <- function(fit, T_shift, meta_info, ndraws = 10) {
  fit_draws <- get_draws_1d_date(fit, "I", ndraws)
  date_mapping <- seq.Date(
    meta_info$T_start_date - T_shift, meta_info$T_end_date,
    by = "1 day"
  )
  fit_draws[, date := date_mapping[as.integer(date)]]
  fit_draws[, c(".chain", ".iteration", "variable") := NULL]
  return(fit_draws[])
}

#' Title
#'
#' @param fit
#' @param data
#' @param meta_info
#' @param ndraws
#'
#' @return
#' @export
#'
#' @examples
summarize_fit <- function(fit, data, meta_info, ndraws = 50) {
  summary <- list()
  T_shift_R <- with(data, L + S + D - G)
  T_shift_latent <- with(data, L + S + D)
  T_shift_onset <- with(data, S)
  T_shift_load <- with(data, 0)

  summary[["R"]] <- get_summary_1d_date(
    fit, "R",
    T_shift = T_shift_R, meta_info = meta_info,
    intervals = c(0.95, 0.5)
  )
  summary[["R"]]$seeding <- FALSE
  summary[["R"]][1:(data$G), "seeding"] <- TRUE

  summary[["R_samples"]] <- get_R_trajectories(
    fit,
    T_shift = T_shift_R, meta_info = meta_info, ndraws = ndraws
  )
  summary[["R_samples"]]$seeding <- FALSE
  summary[["R_samples"]][1:(data$G * ndraws), "seeding"] <- TRUE

  summary[["expected_infections"]] <- get_summary_1d_date(
    fit, "iota",
    T_shift = T_shift_latent, meta_info = meta_info,
    intervals = c(0.95, 0.5)
  )
  summary[["expected_infections"]]$seeding <- FALSE
  summary[["expected_infections"]][1:(data$G * 2), "seeding"] <- TRUE

  if (data$I_sample) {
    summary[["infections"]] <- get_summary_1d_date(
      fit, "I",
      T_shift = T_shift_latent, meta_info = meta_info,
      intervals = c(0.95, 0.5)
    )
    summary[["infections"]]$seeding <- FALSE
    summary[["infections"]][1:(data$G * 2), "seeding"] <- TRUE

    summary[["infections_samples"]] <- get_I_trajectories(
      fit,
      T_shift = T_shift_latent, meta_info = meta_info, ndraws = ndraws
    )
    summary[["infections_samples"]]$seeding <- FALSE
    summary[["infections_samples"]][
      1:(data$G * 2 * ndraws), "seeding"
    ] <- TRUE
  }

  summary[["expected_load"]] <- get_summary_1d_date_log(
    fit, "kappa_log",
    T_shift = T_shift_load, meta_info = meta_info,
    intervals = c(0.95, 0.5)
  )

  summary[["expected_concentration"]] <- get_summary_1d_date_log(
    fit, "pi_log",
    T_shift = T_shift_load, meta_info = meta_info,
    intervals = c(0.95, 0.5)
  )

  summary[["concentration"]] <- get_summary_1d_date(
    fit, "predicted_concentration",
    T_shift = T_shift_load, meta_info = meta_info,
    intervals = c(0.95, 0.5)
  )

  if (data$K > 0) {
    # here we exponentiate to get the multiplicative effect
    summary[["sample_effects"]] <- get_summary_vector_log(
      fit, "eta", colnames(data$X),
      intervals = c(0.95, 0.5)
    )
  }

  return(summary)
}
