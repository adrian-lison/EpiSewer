modeldata_descriptions <- function() {
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

modeldata_var_requirements <- function() {
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

modeldata_defaults <- function() {
  defaults <- list(
    "K" = 0,
    "X" = numeric(0),
    "eta_prior" = numeric(0),
    "init$eta" = numeric(0)
  )
  return(defaults)
}

modeldata_init <- function() {
  modeldata <- list()
  modeldata$init <- list()
  modeldata$meta_info <- list()
  return(modeldata)
}

#' Update meta information in modeldata based on available variables
#'
#' This update function is for all meta information that depends on modeldata
#' from several functions. There are also functions which add meta information
#' directly, in particular when this meta information cannot be solely inferred
#' from modeldata. This function is designed such that calling it will never do
#' harm to the modeldata object and not throw errors if something is missing in
#' the modeldata object.
modeldata_update_metainfo <- function(modeldata) {
  if (modeldata_check(modeldata,
    required = c("L", "S", "T"),
    throw_error = FALSE
  )) {
    modeldata$meta_info$length_I <- with(modeldata, L + S + T)
  }
  if (modeldata_check(modeldata,
    required = c("L", "S", "T", "G"),
    throw_error = FALSE
  )) {
    modeldata$meta_info$length_R <- with(modeldata, L + S + T - G)
  }
  if (modeldata_check(
    modeldata,
    required = c(
      "measured_concentrations",
      "flow",
      "meta_info$composite_window",
      "meta_info$load_per_case"
    ),
    throw_error = FALSE
  )) {
    # crude descriptive estimate of cases at start of time series
    modeldata$meta_info$initial_cases_crude <-
      with(
        modeldata,
        0.1 + measured_concentrations[1] *
          mean(flow[1:meta_info$composite_window]) / meta_info$load_per_case
      )
  }
  return(modeldata)
}

#' Title
#'
#' @param modeldata
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
  function(measurements,
           composite_window = 1,
           date_col = "date",
           measurement_col = "concentration",
           modeldata = modeldata_init()) {
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

    modeldata$T <-
      as.integer(max(measurements[[date_col]]) -
        min(measurements[[date_col]]) + composite_window)
    modeldata$meta_info$T_start_date <-
      min(measurements[[date_col]]) - composite_window + 1
    modeldata$meta_info$T_end_date <- max(measurements[[date_col]])

    modeldata$w <- composite_window
    modeldata$meta_info$composite_window <- composite_window

    measured <- !is.na(measurements[[measurement_col]])
    modeldata$n_measured <- sum(measured)
    modeldata$measured_concentrations <-
      measurements[[measurement_col]][measured]
    modeldata$measured_dates <-
      as.integer(measurements[[date_col]][measured] -
        modeldata$meta_info$T_start_date + 1)
    modeldata$meta_info$measured_dates <-
      measurements[[date_col]][measured]

    return(modeldata)
  }

#' Title
#'
#' @param modeldata
#' @param load_per_case
#'
#' @return
#' @export
#'
#' @examples
load_per_case_assume <-
  function(load_per_case = NULL, modeldata = modeldata_init()) {
    if (is.null(load_per_case)) {
      load_per_case <- tryCatch(
        get_from_env("assumptions", "load_per_case"),
        error = abort_f("Please supply an assumed average shedding load per person.")
      )
    }
    modeldata$meta_info$load_per_case <- load_per_case
    return(modeldata)
  }


#' Title
#'
#' @param modeldata
#' @param flows
#' @param date_col
#' @param flow_col
#'
#' @return
#' @export
#'
#' @examples
flows_observe <-
  function(flows,
           date_col = "date",
           flow_col = "flow",
           modeldata = modeldata_init()) {
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

    modeldata <- tbc(
      "flow_data",
      {
        all_dates <-
          seq.Date(
            modeldata$meta_info$T_start_date,
            modeldata$meta_info$T_end_date,
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
          flows[flows[[date_col]] >= modeldata$meta_info$T_start_date &
            flows[[date_col]] <= modeldata$meta_info$T_end_date, ]
        flows <- flows[order(flows[[date_col]]), ]
        modeldata$flow <- flows[[flow_col]]
      },
      required = c("meta_info$T_start_date", "meta_info$T_end_date")
    )

    return(modeldata)
  }

#' Title
#'
#' @param modeldata
#' @param generation_dist
#'
#' @return
#' @export
#'
#' @examples
generation_dist_assume <-
  function(generation_dist = NULL, modeldata = modeldata_init()) {
    if (is.null(generation_dist)) {
      generation_dist <- tryCatch(
        get_from_env("assumptions", "generation_dist"),
        error = abort_f("Please supply an assumed incubation period distribution.")
      )
    }
    modeldata$G <- length(generation_dist)
    modeldata$generation_dist <- generation_dist
    modeldata$meta_info$length_seeding <- length(generation_dist)
    return(modeldata)
  }

#' Title
#'
#' @param modeldata
#' @param incubation_dist
#'
#' @return
#' @export
#'
#' @examples
incubation_dist_assume <-
  function(incubation_dist = NULL, modeldata = modeldata_init()) {
    if (is.null(incubation_dist)) {
      incubation_dist <- tryCatch(
        get_from_env("assumptions", "incubation_dist"),
        error = abort_f("Please supply an assumed incubation period distribution.")
      )
    }
    modeldata$L <- length(incubation_dist) - 1
    modeldata$incubation_dist <- incubation_dist
    return(modeldata)
  }

#' Title
#'
#' @param modeldata
#' @param shedding_dist
#'
#' @return
#' @export
#'
#' @examples
shedding_dist_assume <-
  function(shedding_dist = NULL, modeldata = modeldata_init()) {
    if (is.null(shedding_dist)) {
      shedding_dist <- tryCatch(
        get_from_env("assumptions", "shedding_dist"),
        error = abort_f("Please supply an assumed shedding load distribution.")
      )
    }
    modeldata$S <- length(shedding_dist) - 1
    # here we account for the scaling factor from cases to load
    modeldata$shedding_dist <-
      tbe(shedding_dist * modeldata$meta_info$load_per_case,
        required = "meta_info$load_per_case"
      )
    return(modeldata)
  }

#' Title
#'
#' @param modeldata
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
    R_level_start_prior =
      set_prior("R_level_start", "normal", mu = 1, sigma = 0.8),
    R_trend_start_prior =
      set_prior("R_trend_start", "normal", mu = 0, sigma = 0.1),
    R_sd_prior =
      set_prior("R_sd", "half-normal", mu = 0, sigma = 0.05),
    ets_diff = FALSE,
    ets_noncentered = TRUE,
    ets_alpha_fixed = NULL,
    ets_alpha_prior = c(50, 50),
    ets_beta_fixed = NULL,
    ets_beta_prior = c(50, 50),
    ets_phi_fixed = 0.9,
    ets_phi_prior = c(50, 5),
    modeldata = modeldata_init()) {

  modeldata$meta_info$R_estimate_approach <- "ets"

  if (is.null(ets_alpha_fixed)) {
    ets_alpha_fixed <- -1
  }
  if (is.null(ets_beta_fixed)) {
    ets_beta_fixed <- -1
  }
  if (is.null(ets_phi_fixed)) {
    ets_phi_fixed <- -1
  }

  modeldata$R_level_start_prior <- R_level_start_prior
  modeldata$R_trend_start_prior <- R_trend_start_prior
  modeldata$R_sd_prior <- R_sd_prior

  modeldata$init$R_level_start <-
    R_level_start_prior$R_level_start_prior[1]
  modeldata$init$R_trend_start <- 1e-4
  modeldata$init$R_sd <- max(R_sd_prior$R_sd_prior[1], 0.1)
  modeldata$init$R_noise <- tbe(
    rep(0, modeldata$meta_info$length_R - 1),
    "meta_info$length_R"
  )

  modeldata$ets_diff <- ets_diff
  modeldata$ets_noncentered <- ets_noncentered

  modeldata$ets_alpha_fixed <- ets_alpha_fixed
  if (ets_alpha_fixed >= 0) {
    modeldata$ets_alpha_prior <- numeric(0)
    modeldata$init$ets_alpha <- numeric(0)
  } else {
    modeldata$ets_alpha_prior <- ets_alpha_prior
    modeldata$init$ets_alpha <- 0.5
  }


  modeldata$ets_beta_fixed <- ets_beta_fixed
  if (ets_beta_fixed >= 0) {
    modeldata$ets_beta_prior <- numeric(0)
    modeldata$init$ets_beta <- numeric(0)
  } else {
    modeldata$ets_beta_prior <- ets_beta_prior
    modeldata$init$ets_beta <- 0.5
  }


  modeldata$ets_phi_fixed <- ets_phi_fixed
  if (ets_phi_fixed >= 0) {
    modeldata$ets_phi_prior <- numeric(0)
    modeldata$init$ets_phi <- numeric(0)
  } else {
    modeldata$ets_phi_prior <- ets_phi_prior
    modeldata$init$ets_phi <- 0.9
  }

  return(modeldata)
}

R_estimate_rw <- function(
    R_start_prior =
      set_prior("R_level_start", "normal", mu = 1, sigma = 0.8),
    R_sd_prior =
      set_prior("R_sd", "half-normal", mu = 0, sigma = 0.05),
    rw_diff = FALSE,
    rw_noncentered = TRUE,
    modeldata = modeldata_init()) {

  modeldata <- R_estimate_ets(
    modeldata = modeldata,
    R_level_start_prior = R_start_prior,
    R_sd_prior = R_sd_prior,
    ets_diff = rw_diff,
    ets_noncentered = rw_noncentered,
    ets_alpha_fixed = 1,
    ets_beta_fixed = 0,
    ets_phi_fixed = 0
    )
  return(modeldata)
}

#' Title
#'
#' @param modeldata
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
    knot_distance = 1,
    spline_degree = 3,
    bs_coeff_ar_start_prior = set_prior("bs_coeff_ar_start",
      "normal",
      mu = 0,
      sigma = 0.5
    ),
    bs_coeff_ar_sd_prior = set_prior("bs_coeff_ar_sd",
      "half-normal",
      mu = 0,
      sigma = 0.2
    ),
    modeldata = modeldata_init()) {
  modeldata$meta_info$R_estimate_approach <- "splines"

  modeldata <- tbc("spline_definition",
    {
      knots <- seq(1, modeldata$meta_info$length_R, by = knot_distance)
      B <-
        splines::bs(
          1:modeldata$meta_info$length_R,
          knots = knots,
          degree = spline_degree,
          intercept = FALSE
        )
      modeldata$meta_info$R_knots <- knots
      modeldata$meta_info$B <- B
      modeldata$bs_n_basis <- ncol(B)
      B_sparse <- suppressMessages(rstan::extract_sparse_parts(B))
      modeldata$bs_n_w <- length(B_sparse$w)
      modeldata$bs_w <- B_sparse$w
      modeldata$bs_v <- B_sparse$v
      modeldata$bs_u <- B_sparse$u
    },
    required = "meta_info$length_R"
  )

  modeldata$bs_coeff_ar_start_prior <- bs_coeff_ar_start_prior
  modeldata$bs_coeff_ar_sd_prior <- bs_coeff_ar_sd_prior

  modeldata$init$bs_coeff_ar_start <- 0
  modeldata$init$bs_coeff_ar_sd <- 0.1
  modeldata$init$bs_coeff_noise <- tbe(
    rep(0, modeldata$meta_info$length_R - 1),
    "meta_info$length_R"
  )

  return(modeldata)
}

#' Title
#'
#' @param modeldata
#' @param iota_log_ar_start_prior
#' @param iota_log_ar_sd_prior
#'
#' @return
#' @export
#'
#' @examples
seeding_estimate <- function(
    iota_log_ar_start_prior = NULL,
    iota_log_ar_sd_prior = set_prior("iota_log_ar_sd",
      "half-normal",
      mu = 0.05,
      sigma = 0.025
    ),
    modeldata = modeldata_init()) {
  new_modeldata <- as.list(environment())
  modeldata <-
    c(modeldata, new_modeldata[names(new_modeldata) != "modeldata"])

  if (is.null(modeldata$iota_log_ar_start_prior)) {
    modeldata$iota_log_ar_start_prior <- tbe(
      set_prior(
        "iota_log_ar_start",
        "normal (mu based on crude empirical estimate of cases)",
        mu = log(modeldata$meta_info$initial_cases_crude),
        sigma = 1
      ),
      "meta_info$initial_cases_crude"
    )
  }

  modeldata$init$iota_log_ar_start <- tbe(
    log(modeldata$meta_info$initial_cases_crude),
    "meta_info$initial_cases_crude"
  )
  modeldata$init$iota_log_ar_sd <- 1
  modeldata$init$iota_log_ar_noise <- tbe(
    rep(0, modeldata$meta_info$length_seeding - 1),
    "meta_info$length_seeding"
  )

  return(modeldata)
}

#' Title
#'
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
infection_noise_none <- function(modeldata = modeldata_init()) {
  modeldata$I_sample <- FALSE
  modeldata$I_overdispersion <- FALSE
  modeldata$init$I <- numeric(0)
  modeldata$init$I_log <- numeric(0)
  return(modeldata)
}

#' Title
#'
#' @param modeldata
#' @param I_overdispersion
#' @param I_xi_prior
#'
#' @return
#' @export
#'
#' @examples
infection_noise_estimate <-
  function(I_overdispersion = FALSE,
           I_xi_prior = set_prior("I_xi", "normal", mu = 0, sigma = 1),
           modeldata = modeldata_init()) {
    new_modeldata <- as.list(environment())
    modeldata <-
      c(modeldata, new_modeldata[names(new_modeldata) != "modeldata"])
    modeldata$I_sample <- TRUE
    modeldata$init$I <- tbe(
      rep(
        modeldata$meta_info$initial_cases_crude,
        modeldata$meta_info$length_I
      ),
      c("meta_info$initial_cases_crude", "meta_info$length_I")
    )
    modeldata$init$I_log <- tbe(
      rep(
        log(modeldata$meta_info$initial_cases_crude),
        modeldata$meta_info$length_I
      ),
      c("meta_info$initial_cases_crude", "meta_info$length_I")
    )

    if (modeldata$I_overdispersion) {
      modeldata$init$I_xi <- 0.05
    } else {
      modeldata$I_xi_prior <- numeric(0)
      modeldata$init$I_xi <- numeric(0)
    }

    return(modeldata)
  }

#' Title
#'
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
sample_effects_none <- function(modeldata = modeldata_init()) {
  modeldata$K <- 0
  modeldata$X <- numeric(0)
  modeldata$eta_prior <- numeric(0)
  modeldata$init$eta <- numeric(0)
  return(modeldata)
}

#' Title
#'
#' @param modeldata
#' @param design_matrix
#' @param eta_prior
#'
#' @return
#' @export
#'
#' @examples
sample_effects_estimate_matrix <-
  function(design_matrix,
           eta_prior = set_prior("eta", "normal", mu = 0, sigma = 1),
           modeldata = modeldata_init()) {
    modeldata <- tbc(
      "check_design_matrix",
      {
        if (!(modeldata$T == nrow(design_matrix))) {
          abort(
            paste(
              "Mismatch: Modeled time period has",
              modeldata$T,
              "days, design matrix for sample date effects has",
              nrow(design_matrix),
              "rows."
            )
          )
        }
      },
      required = "T"
    )
    modeldata$K <- ncol(design_matrix)
    modeldata$X <- design_matrix

    modeldata$eta_prior <- eta_prior

    modeldata$init$eta <- rep(0, modeldata$K)

    return(modeldata)
  }

#' Title
#'
#' @param modeldata
#' @param eta_prior
#'
#' @return
#' @export
#'
#' @examples
sample_effects_estimate_weekday <-
  function(eta_prior = set_prior("eta", "normal", mu = 0, sigma = 1),
           modeldata = modeldata_init()) {
    modeldata <- tbc(
      "weekday_design_matrix",
      {
        weekdays <- lubridate::wday(
          seq.Date(
            modeldata$meta_info$T_start_date,
            modeldata$meta_info$T_end_date,
            by = "1 day"
          ),
          label = TRUE
        )
        design_matrix <- model.matrix(
          ~wday,
          data.frame(wday = weekdays),
          contrasts.arg = list(wday = "contr.treatment")
        )[, -1]
        modeldata <-
          sample_effects_estimate_matrix(modeldata, design_matrix, eta_prior)
      },
      required = c("meta_info$T_start_date", "meta_info$T_end_date")
    )
  }

#' Title
#'
#' @param modeldata
#' @param sigma_prior
#'
#' @return
#' @export
#'
#' @examples
measurement_noise_estimate <-
  function(sigma_prior = set_prior("sigma", "normal", mu = 0, sigma = 1),
           modeldata = modeldata_init()) {
    new_modeldata <- as.list(environment())
    modeldata <-
      c(modeldata, new_modeldata[names(new_modeldata) != "modeldata"])

    modeldata$init$sigma <-
      0.1 # roughly corresponds to a 10% coefficient of variation

    return(modeldata)
  }
