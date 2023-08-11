standata_descriptions <- function() {
  descriptions <- list(
    "meta_info$composite_window" =
      "window length for composite samples in days",
    "meta_info$length_seeding" =
      "length of seeding phase for infections",
    "meta_info$length_I" =
      "number of days over which infections are modeled",
    "meta_info$length_R" =
      "number of days over which Rt is modeled",
    "meta_info$load_per_case" =
      "assumed overall load shed per individual",
    "meta_info$initial_cases_crude" =
      "empirical estimate for #cases at start of time period"
  )
  return(descriptions)
}

standata_var_requirements <- function() {
  requirements <- list(
    "meta_info$initial_cases_crude" = "flows_observe",
    "meta_info$length_seeding" = "generation_dist_assume",
    "meta_info$length_I" = c("incubation_dist_assume", "shedding_dist_assume"),
    "meta_info$length_R" = c(
      "incubation_dist_assume",
      "shedding_dist_assume",
      "generation_dist_assume"
    )
  )
  return(requirements)
}

standata_defaults <- function() {
  defaults <- list(
    "K" = 0,
    "X" = numeric(0),
    "eta_prior" = numeric(0),
    "init$eta" = numeric(0)
  )
  return(defaults)
}

standata_init <- function() {
  standata <- list()
  standata$init <- list()
  standata$meta_info <- list()
  return(standata)
}

#' Update meta information in standata based on available variables
#'
#' This update function is for all meta information that depends on standata
#' from several functions. There are also functions which add meta information
#' directly, in particular when this meta information cannot be solely inferred
#' from standata. This function is designed such that calling it will never do
#' harm to the standata object and not throw errors if something is missing in
#' the standata object.
standata_update_metainfo <- function(standata) {
  if (standata_check(standata,
    required = c("L", "S", "T"),
    throw_error = F
  )) {
    standata$meta_info$length_I <- with(standata, L + S + T)
  }
  if (standata_check(standata,
    required = c("L", "S", "T", "G"),
    throw_error = F
  )) {
    standata$meta_info$length_R <- with(standata, L + S + T - G)
  }
  if (standata_check(
    standata,
    required = c(
      "measured_concentrations",
      "flow",
      "meta_info$composite_window",
      "meta_info$load_per_case"
    ),
    throw_error = F
  )) {
    # crude descriptive estimate of cases at start of time series
    standata$meta_info$initial_cases_crude <-
      with(
        standata,
        0.1 + measured_concentrations[1] *
          mean(flow[1:meta_info$composite_window]) / meta_info$load_per_case
      )
  }
  return(standata)
}

#' Title
#'
#' @param standata
#' @param measurements
#' @param composite_window
#' @param date_col
#' @param measurement_col
#'
#' @return
#' @export
#'
#' @examples
measurements_observe <-
  function(standata = standata_init(),
           measurements,
           composite_window = 1,
           date_col = "date",
           measurement_col = "concentration") {
    if (is.null(measurements)) {
      abort("Please supply measurement data.")
    }

    required_data_cols <- c(date_col, measurement_col)
    if (!all(required_data_cols %in% names(measurements))) {
      abort(
        paste(
          "The following columns must be present",
          "in the provided measurements `data.frame`:",
          paste(required_data_cols, collapse = ", ")
        )
      )
    }

    measurements[[date_col]] <- as.Date(measurements[[date_col]])

    standata$T <-
      as.integer(max(measurements[[date_col]]) -
        min(measurements[[date_col]]) + composite_window)
    standata$meta_info$T_start_date <-
      min(measurements[[date_col]]) - composite_window + 1
    standata$meta_info$T_end_date <- max(measurements[[date_col]])

    standata$w <- composite_window
    standata$meta_info$composite_window <- composite_window

    measured <- !is.na(measurements[[measurement_col]])
    standata$n_measured <- sum(measured)
    standata$measured_concentrations <-
      measurements[[measurement_col]][measured]
    standata$measured_dates <-
      as.integer(measurements[[date_col]][measured] -
        standata$meta_info$T_start_date + 1)
    standata$meta_info$measured_dates <-
      measurements[[date_col]][measured]

    return(standata)
  }

#' Title
#'
#' @param standata
#' @param load_per_case
#'
#' @return
#' @export
#'
#' @examples
load_per_case_assume <-
  function(standata = standata_init(), load_per_case) {
    if (is.null(load_per_case)) {
      abort("Please supply an assumed average shedding load per person.")
    }
    standata$meta_info$load_per_case <- load_per_case
    return(standata)
  }


#' Title
#'
#' @param standata
#' @param flows
#' @param date_col
#' @param flow_col
#'
#' @return
#' @export
#'
#' @examples
flows_observe <-
  function(standata = standata_init(),
           flows,
           date_col = "date",
           flow_col = "flow") {
    if (is.null(flows)) {
      abort("Please supply flow data.")
    }

    required_data_cols <- c(date_col, flow_col)
    if (!all(required_data_cols %in% names(flows))) {
      abort(
        paste(
          "The following columns must be present",
          "in the provided flow `data.frame`:",
          paste(required_data_cols, collapse = ", ")
        )
      )
    }

    flows[[date_col]] <- as.Date(flows[[date_col]])

    if (any(duplicated(flows[[date_col]]))) {
      abort("Flow data is ambigious, duplicate dates found.")
    }

    standata <- tbef(
      "flow_data",
      {
        all_dates <-
          seq.Date(
            standata$meta_info$T_start_date,
            standata$meta_info$T_end_date,
            by = "1 day"
          )
        missing_flow_dates <-
          as.Date(
            setdiff(all_dates, flows[[date_col]][!is.na(flows[[flow_col]])]),
            origin = lubridate::origin
          )
        if (length(missing_flow_dates) > 0) {
          abort(paste(
            "Missing flow values for the following dates:",
            paste(missing_flow_dates, collapse = ", ")
          ))
        }
        flows <-
          flows[flows[[date_col]] >= standata$meta_info$T_start_date &
            flows[[date_col]] <= standata$meta_info$T_end_date, ]
        flows <- flows[order(flows[[date_col]]), ]
        standata$flow <- flows[[flow_col]]
      },
      required = c("meta_info$T_start_date", "meta_info$T_end_date")
    )

    return(standata)
  }

#' Title
#'
#' @param standata
#' @param generation_dist
#'
#' @return
#' @export
#'
#' @examples
generation_dist_assume <-
  function(standata = standata_init(), generation_dist) {
    if (is.null(generation_dist)) {
      abort("Please supply an assumed generation time distribution.")
    }
    standata$G <- length(generation_dist)
    standata$generation_dist <- generation_dist
    standata$meta_info$length_seeding <- length(generation_dist)
    return(standata)
  }

#' Title
#'
#' @param standata
#' @param incubation_dist
#'
#' @return
#' @export
#'
#' @examples
incubation_dist_assume <-
  function(standata = standata_init(), incubation_dist) {
    if (is.null(incubation_dist)) {
      abort("Please supply an assumed incubation period distribution.")
    }
    standata$L <- length(incubation_dist) - 1
    standata$incubation_dist <- incubation_dist
    return(standata)
  }

#' Title
#'
#' @param standata
#' @param shedding_dist
#'
#' @return
#' @export
#'
#' @examples
shedding_dist_assume <-
  function(standata = standata_init(), shedding_dist) {
    if (is.null(shedding_dist)) {
      abort("Please supply an assumed shedding load distribution.")
    }
    standata$S <- length(shedding_dist) - 1
    # here we account for the scaling factor from cases to load
    standata$shedding_dist <-
      tbe(shedding_dist * standata$meta_info$load_per_case,
        required = "meta_info$load_per_case"
      )
    return(standata)
  }

#' Title
#'
#' @param standata
#' @param R_level_start_prior
#' @param R_trend_start_prior
#' @param R_sd_prior
#' @param ets_diff
#' @param ets_noncentered
#' @param ets_alpha_fixed
#' @param ets_alpha_prior
#' @param ets_beta_fixed
#' @param ets_beta_prior
#' @param ets_phi_fixed
#' @param ets_phi_prior
#'
#' @return
#' @export
#'
#' @examples
R_estimate_ets <- function(
    standata = standata_init(),
    R_level_start_prior =
      stan_prior("R_level_start", "normal", mu = 1, sigma = 0.8),
    R_trend_start_prior =
      stan_prior("R_trend_start", "normal", mu = 0, sigma = 0.1),
    R_sd_prior =
      stan_prior("R_sd", "half-normal", mu = 0, sigma = 0.05),
    ets_diff = FALSE,
    ets_noncentered = TRUE,
    ets_alpha_fixed = NULL,
    ets_alpha_prior = c(50, 50),
    ets_beta_fixed = NULL,
    ets_beta_prior = c(50, 50),
    ets_phi_fixed = 0.9,
    ets_phi_prior = c(50, 5)) {

  standata$meta_info$R_estimate_approach <- "ets"

  if (is.null(ets_alpha_fixed)) {
    ets_alpha_fixed <- -1
  }
  if (is.null(ets_beta_fixed)) {
    ets_beta_fixed <- -1
  }
  if (is.null(ets_phi_fixed)) {
    ets_phi_fixed <- -1
  }

  standata$R_level_start_prior <- R_level_start_prior
  standata$R_trend_start_prior <- R_trend_start_prior
  standata$R_sd_prior <- R_sd_prior

  standata$init$R_level_start <-
    R_level_start_prior$R_level_start_prior[1]
  standata$init$R_trend_start <- 1e-4
  standata$init$R_sd <- max(R_sd_prior$R_sd_prior[1], 0.1)
  standata$init$R_noise <- tbe(
    rep(0, standata$meta_info$length_R - 1),
    "meta_info$length_R"
  )

  standata$ets_diff <- ets_diff
  standata$ets_noncentered <- ets_noncentered

  standata$ets_alpha_fixed <- ets_alpha_fixed
  if (ets_alpha_fixed >= 0) {
    standata$ets_alpha_prior <- numeric(0)
  } else {
    standata$ets_alpha_prior <- ets_alpha_prior
  }
  standata$init$ets_alpha <- 0.5

  standata$ets_beta_fixed <- ets_beta_fixed
  if (ets_beta_fixed >= 0) {
    standata$ets_beta_prior <- numeric(0)
  } else {
    standata$ets_beta_prior <- ets_beta_prior
  }
  standata$init$ets_beta <- 0.5

  standata$ets_phi_fixed <- ets_phi_fixed
  if (ets_phi_fixed >= 0) {
    standata$ets_phi_prior <- numeric(0)
  } else {
    standata$ets_phi_prior <- ets_phi_prior
  }
  standata$init$ets_phi <- 0.9

  return(standata)
}

R_estimate_rw <- function(
    standata = standata_init(),
    R_start_prior =
      stan_prior("R_level_start", "normal", mu = 1, sigma = 0.8),
    R_sd_prior =
      stan_prior("R_sd", "half-normal", mu = 0, sigma = 0.05),
    rw_diff = FALSE,
    rw_noncentered = TRUE) {

  standata <- R_estimate_ets(
    standata = standata,
    R_level_start_prior = R_start_prior,
    R_sd_prior = R_sd_prior,
    ets_diff = rw_diff,
    ets_noncentered = rw_noncentered,
    ets_alpha_fixed = 1,
    ets_beta_fixed = 0,
    ets_phi_fixed = 0
    )
  return(standata)
}

#' Title
#'
#' @param standata
#' @param knot_distance
#' @param spline_degree
#' @param bs_coeff_ar_start_prior
#' @param bs_coeff_ar_sd_prior
#'
#' @return
#' @export
#'
#' @examples
R_estimate_splines <- function(
    standata = standata_init(),
    knot_distance = 1,
    spline_degree = 3,
    bs_coeff_ar_start_prior = stan_prior("bs_coeff_ar_start",
      "normal",
      mu = 0,
      sigma = 0.5
    ),
    bs_coeff_ar_sd_prior = stan_prior("bs_coeff_ar_sd",
      "half-normal",
      mu = 0,
      sigma = 0.2
    )) {
  standata$meta_info$R_estimate_approach <- "splines"

  standata <- tbef("spline_definition",
    {
      knots <- seq(1, standata$meta_info$length_R, by = knot_distance)
      B <-
        splines::bs(
          1:standata$meta_info$length_R,
          knots = knots,
          degree = spline_degree,
          intercept = F
        )
      standata$meta_info$R_knots <- knots
      standata$meta_info$B <- B
      standata$bs_n_basis <- ncol(B)
      B_sparse <- suppressMessages(rstan::extract_sparse_parts(B))
      standata$bs_n_w <- length(B_sparse$w)
      standata$bs_w <- B_sparse$w
      standata$bs_v <- B_sparse$v
      standata$bs_u <- B_sparse$u
    },
    required = "meta_info$length_R"
  )

  standata$bs_coeff_ar_start_prior <- bs_coeff_ar_start_prior
  standata$bs_coeff_ar_sd_prior <- bs_coeff_ar_sd_prior

  standata$init$bs_coeff_ar_start <- 0
  standata$init$bs_coeff_ar_sd <- 0.1
  standata$init$bs_coeff_noise <- tbe(
    rep(0, standata$meta_info$length_R - 1),
    "meta_info$length_R"
  )

  return(standata)
}

#' Title
#'
#' @param standata
#' @param iota_log_ar_start_prior
#' @param iota_log_ar_sd_prior
#'
#' @return
#' @export
#'
#' @examples
seeding_estimate <- function(
    standata = standata_init(),
    iota_log_ar_start_prior = NULL,
    iota_log_ar_sd_prior = stan_prior("iota_log_ar_sd",
      "half-normal",
      mu = 0.05,
      sigma = 0.025
    )) {
  new_standata <- as.list(environment())
  standata <-
    c(standata, new_standata[names(new_standata) != "standata"])

  if (is.null(standata$iota_log_ar_start_prior)) {
    standata$iota_log_ar_start_prior <- tbe(
      stan_prior(
        "iota_log_ar_start",
        "normal (mu based on crude empirical estimate of cases)",
        mu = log(standata$meta_info$initial_cases_crude),
        sigma = 1
      ),
      "meta_info$initial_cases_crude"
    )
  }

  standata$init$iota_log_ar_start <- tbe(
    log(standata$meta_info$initial_cases_crude),
    "meta_info$initial_cases_crude"
  )
  standata$init$iota_log_ar_sd <- 1
  standata$init$iota_log_ar_noise <- tbe(
    rep(0, standata$meta_info$length_seeding - 1),
    "meta_info$length_seeding"
  )

  return(standata)
}

#' Title
#'
#' @param standata
#'
#' @return
#' @export
#'
#' @examples
infection_noise_none <- function(standata = standata_init()) {
  standata$I_sample <- FALSE
  standata$I_overdispersion <- FALSE
  standata$init$I <- numeric(0)
  standata$init$I_log <- numeric(0)
  return(standata)
}

#' Title
#'
#' @param standata
#' @param I_overdispersion
#' @param I_xi_prior
#'
#' @return
#' @export
#'
#' @examples
infection_noise_estimate <-
  function(standata = standata_init(),
           I_overdispersion = FALSE,
           I_xi_prior = stan_prior("I_xi", "normal", mu = 0, sigma = 1)) {
    new_standata <- as.list(environment())
    standata <-
      c(standata, new_standata[names(new_standata) != "standata"])
    standata$I_sample <- TRUE
    standata$init$I <- tbe(
      rep(
        standata$meta_info$initial_cases_crude,
        standata$meta_info$length_I
      ),
      c("meta_info$initial_cases_crude", "meta_info$length_I")
    )
    standata$init$I_log <- tbe(
      rep(
        log(standata$meta_info$initial_cases_crude),
        standata$meta_info$length_I
      ),
      c("meta_info$initial_cases_crude", "meta_info$length_I")
    )

    if (standata$I_overdispersion) {
      standata$init$I_xi <- 0.05
    } else {
      standata$I_xi_prior <- numeric(0)
      standata$init$I_xi <- numeric(0)
    }

    return(standata)
  }

#' Title
#'
#' @param standata
#'
#' @return
#' @export
#'
#' @examples
sample_effects_none <- function(standata = standata_init()) {
  standata$K <- 0
  standata$X <- numeric(0)
  standata$eta_prior <- numeric(0)
  standata$init$eta <- numeric(0)
  return(standata)
}

#' Title
#'
#' @param standata
#' @param design_matrix
#' @param eta_prior
#'
#' @return
#' @export
#'
#' @examples
sample_effects_estimate_matrix <-
  function(standata = standata_init(),
           design_matrix,
           eta_prior = stan_prior("eta", "normal", mu = 0, sigma = 1)) {
    standata <- tbef(
      "check_design_matrix",
      {
        if (!(standata$T == nrow(design_matrix))) {
          abort(
            paste(
              "Mismatch: Modeled time period has",
              standata$T,
              "days, design matrix for sample date effects has",
              nrow(design_matrix),
              "rows."
            )
          )
        }
      },
      required = "T"
    )
    standata$K <- ncol(design_matrix)
    standata$X <- design_matrix

    standata$eta_prior <- eta_prior

    standata$init$eta <- rep(0, standata$K)

    return(standata)
  }

#' Title
#'
#' @param standata
#' @param eta_prior
#'
#' @return
#' @export
#'
#' @examples
sample_effects_estimate_weekday <-
  function(standata = standata_init(),
           eta_prior = stan_prior("eta", "normal", mu = 0, sigma = 1)) {
    standata <- tbef(
      "weekday_design_matrix",
      {
        weekdays <- lubridate::wday(
          seq.Date(
            standata$meta_info$T_start_date,
            standata$meta_info$T_end_date,
            by = "1 day"
          ),
          label = T
        )
        design_matrix <- model.matrix(
          ~wday,
          data.frame(wday = weekdays),
          contrasts.arg = list(wday = "contr.treatment")
        )[, -1]
        standata <-
          sample_effects_estimate_matrix(standata, design_matrix, eta_prior)
      },
      required = c("meta_info$T_start_date", "meta_info$T_end_date")
    )
  }

#' Title
#'
#' @param standata
#' @param sigma_prior
#'
#' @return
#' @export
#'
#' @examples
measurement_noise_estimate <-
  function(standata = standata_init(),
           sigma_prior = stan_prior("sigma", "normal", mu = 0, sigma = 1)) {
    new_standata <- as.list(environment())
    standata <-
      c(standata, new_standata[names(new_standata) != "standata"])

    standata$init$sigma <-
      0.1 # roughly corresponds to a 10% coefficient of variation

    return(standata)
  }
