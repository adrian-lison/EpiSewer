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
  modeldata$checks <- list()
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
#' @param data
#' @param composite_window
#' @param date_col
#' @param concentration_col
#' @param replicate_col
#'
#' @return
#' @export
#'
#' @examples
concentrations_observe <-
  function(data,
           composite_window = 1,
           date_col = "date",
           concentration_col = "concentration",
           replicate_col = NULL,
           modeldata = modeldata_init()) {
    if (is.null(data)) {
      abort("Please supply measurement data.")
    }

    required_data_cols <- c(date_col, concentration_col, replicate_col)
    if (!all(required_data_cols %in% names(data))) {
      abort(
        paste(
          "The following columns must be present",
          "in the provided measurements `data.frame`:",
          paste(required_data_cols, collapse = ", ")
        )
      )
    }

    data[[date_col]] <- as.Date(data[[date_col]])

    if (is.null(replicate_col)) {
      if (any(duplicated(data[[date_col]]))) {
        abort(
          paste(
            "Duplicated dates found in measurements `data.frame`.",
            "If your data contains replicate measurements, please add a column",
            "with replicate IDs to your `data.frame` and provide its name via",
            "the `replicate_col` argument."
          )
        )
      }
    }

    modeldata$T <-
      as.integer(max(data[[date_col]]) -
        min(data[[date_col]]) + composite_window)

    modeldata$meta_info$T_start_date <-
      min(data[[date_col]]) - composite_window + 1
    modeldata$meta_info$T_end_date <- max(data[[date_col]])

    modeldata$w <- composite_window
    modeldata$meta_info$composite_window <- composite_window


    measured <- !is.na(data[[concentration_col]])
    modeldata$n_measured <- sum(measured)
    modeldata$n_samples <- length(unique(data[[date_col]][measured]))
    measured_dates <- as.integer(
      data[[date_col]][measured] - modeldata$meta_info$T_start_date + 1
    )
    modeldata$sample_to_date <- sort(unique(measured_dates))
    modeldata$measure_to_sample <- sapply(
      measured_dates,
      function(x) which(x == modeldata$sample_to_date)[[1]]
    )
    modeldata$measured_concentrations <- data[[concentration_col]][measured]
    modeldata$meta_info$measured_dates <- data[[date_col]][measured]

    if (!is.null(replicate_col)) {
      modeldata$replicate_ids <- as.integer(data[[replicate_col]][measured])
    }

    return(modeldata)
  }

droplets_observe <-
  function(data,
           composite_window = 1,
           date_col = "date",
           droplets_col = "droplets",
           replicate_col = NULL,
           modeldata = modeldata_init()) {
    abort(paste(
      "Specification of measurements via ddPCR droplet count",
      "is not implemented yet."
    ))
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
#' @param data
#' @param date_col
#' @param flow_col
#'
#' @return
#' @export
#'
#' @examples
flows_observe <-
  function(data,
           date_col = "date",
           flow_col = "flow",
           modeldata = modeldata_init()) {
    if (is.null(data)) {
      abort("Please supply flow data.")
    }

    required_data_cols <- c(date_col, flow_col)
    if (!all(required_data_cols %in% names(data))) {
      abort(
        paste(
          "The following columns must be present",
          "in the provided flow `data.frame`:",
          paste(required_data_cols, collapse = ", ")
        )
      )
    }

    data[[date_col]] <- as.Date(data[[date_col]])

    if (any(duplicated(data[[date_col]]))) {
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
            setdiff(all_dates, data[[date_col]][!is.na(data[[flow_col]])]),
            origin = lubridate::origin
          )
        if (length(missing_flow_dates) > 0) {
          abort(paste(
            "Missing flow values for the following dates:",
            paste(missing_flow_dates, collapse = ", ")
          ))
        }
        data <-
          data[data[[date_col]] >= modeldata$meta_info$T_start_date &
            data[[date_col]] <= modeldata$meta_info$T_end_date, ]
        data <- data[order(data[[date_col]]), ]
        modeldata$flow <- data[[flow_col]]
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
#' @param level_prior_mu
#' @param level_prior_sigma
#' @param trend_prior_mu
#' @param trend_prior_sigma
#' @param sd_prior_mu
#' @param sd_prior_sigma
#' @param smooth_fixed
#' @param smooth_prior_shapes
#' @param trend_smooth_fixed
#' @param trend_smooth_prior_shapes
#' @param dampen_prior_shapes
#' @param dampen_fixed
#' @param differenced
#' @param noncentered
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
R_estimate_ets <- function(
    level_prior_mu = 1,
    level_prior_sigma = 0.8,
    trend_prior_mu = 0,
    trend_prior_sigma = 0.1,
    sd_prior_mu = 0,
    sd_prior_sigma = 0.05,
    smooth_fixed = NULL,
    smooth_prior_shapes = c(50, 50),
    trend_smooth_fixed = NULL,
    trend_smooth_prior_shapes = c(50, 50),
    dampen_fixed = 0.9,
    dampen_prior_shapes = c(50, 5),
    differenced = FALSE,
    noncentered = TRUE,
    modeldata = modeldata_init()) {
  modeldata$meta_info$R_estimate_approach <- "ets"

  modeldata$R_level_start_prior <- set_prior(
    "R_level_start", "normal",
    mu = level_prior_mu, sigma = level_prior_sigma
  )

  modeldata$R_trend_start_prior <- set_prior(
    "R_trend_start", "normal",
    mu = trend_prior_mu, sigma = trend_prior_sigma
  )

  modeldata$R_sd_prior <- set_prior(
    "R_sd", "truncated normal",
    mu = sd_prior_mu, sigma = sd_prior_sigma
  )

  modeldata$init$R_level_start <-
    modeldata$R_level_start_prior$R_level_start_prior[1]
  modeldata$init$R_trend_start <- 1e-4
  modeldata$init$R_sd <- max(modeldata$R_sd_prior$R_sd_prior[1], 0.1)
  modeldata$init$R_noise <- tbe(
    rep(0, modeldata$meta_info$length_R - 1),
    "meta_info$length_R"
  )

  if (is.null(smooth_fixed)) {
    smooth_fixed <- -1
  }
  modeldata$ets_alpha_fixed <- smooth_fixed
  if (smooth_fixed >= 0) {
    modeldata$ets_alpha_prior <- numeric(0)
    modeldata$init$ets_alpha <- numeric(0)
  } else {
    modeldata$ets_alpha_prior <- set_prior("ets_alpha", "beta",
      alpha = smooth_prior_shapes[1], beta = smooth_prior_shapes[2]
    )
    modeldata$init$ets_alpha <- 0.5
  }

  if (is.null(trend_smooth_fixed)) {
    trend_smooth_fixed <- -1
  }
  modeldata$ets_beta_fixed <- trend_smooth_fixed
  if (trend_smooth_fixed >= 0) {
    modeldata$ets_beta_prior <- numeric(0)
    modeldata$init$ets_beta <- numeric(0)
  } else {
    modeldata$ets_beta_prior <- set_prior("ets_beta", "beta",
      alpha = trend_smooth_prior_shapes[1], beta = trend_smooth_prior_shapes[2]
    )
    modeldata$init$ets_beta <- 0.5
  }

  if (is.null(dampen_fixed)) {
    dampen_fixed <- -1
  }
  modeldata$ets_phi_fixed <- dampen_fixed
  if (dampen_fixed >= 0) {
    modeldata$ets_phi_prior <- numeric(0)
    modeldata$init$ets_phi <- numeric(0)
  } else {
    modeldata$ets_phi_prior <- set_prior("ets_phi", "beta",
      alpha = dampen_prior_shapes[1], beta = dampen_prior_shapes[2]
    )
    modeldata$init$ets_phi <- 0.9
  }

  modeldata$ets_diff <- differenced
  modeldata$ets_noncentered <- noncentered

  return(modeldata)
}

#' Title
#'
#' @param level_prior_mu
#' @param level_prior_sigma
#' @param sd_prior_mu
#' @param sd_prior_sigma
#' @param differenced
#' @param noncentered
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
R_estimate_rw <- function(
    intercept_prior_mu = 1,
    intercept_prior_sigma = 0.8,
    sd_prior_mu = 0,
    sd_prior_sigma = 0.05,
    differenced = FALSE,
    noncentered = TRUE,
    modeldata = modeldata_init()) {
  modeldata <- R_estimate_ets(
    level_prior_mu = intercept_prior_mu,
    level_prior_sigma = intercept_prior_sigma,
    sd_prior_mu = sd_prior_mu,
    sd_prior_sigma = sd_prior_sigma,
    smooth_fixed = 1,
    trend_smooth_fixed = 0,
    dampen_fixed = 0,
    differenced = differenced,
    noncentered = noncentered,
    modeldata = modeldata
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
    coef_intercept_prior_mu = 0,
    coef_intercept_prior_sigma = 0.5,
    coef_sd_prior_mu = 0,
    coef_sd_prior_sigma = 0.2,
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

      modeldata$init$bs_coeff_ar_start <- 0
      modeldata$init$bs_coeff_ar_sd <- 0.1
      modeldata$init$bs_coeff_noise <- rep(0, modeldata$bs_n_basis - 1)
    },
    required = "meta_info$length_R"
  )

  modeldata$bs_coeff_ar_start_prior <- set_prior("bs_coeff_ar_start",
    "normal",
    mu = coef_intercept_prior_mu,
    sigma = coef_intercept_prior_sigma
  )
  modeldata$bs_coeff_ar_sd_prior <- set_prior("bs_coeff_ar_sd",
    "truncated normal",
    mu = coef_sd_prior_mu,
    sigma = coef_sd_prior_sigma
  )

  return(modeldata)
}

#' Title
#'
#' @param intercept_prior_mu
#' @param intercept_prior_sigma
#' @param sd_prior_mu
#' @param sd_prior_sigma
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
seeding_estimate <- function(
    intercept_prior_mu = NULL,
    intercept_prior_sigma = NULL,
    sd_prior_mu = 0.05,
    sd_prior_sigma = 0.025,
    modeldata = modeldata_init()) {
  if (!(is.null(intercept_prior_mu) || is.null(intercept_prior_sigma))) {
    modeldata$iota_log_ar_start_prior <- set_prior(
      "iota_log_ar_start",
      "normal",
      mu = intercept_prior_mu,
      sigma = intercept_prior_sigma
    )
  }

  modeldata$iota_log_ar_sd_prior <- set_prior(
    "iota_log_ar_sd",
    "truncated normal",
    mu = sd_prior_mu,
    sigma = sd_prior_sigma
  )

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
#' @param overdispersion
#' @param overdispersion_prior_mu
#' @param overdispersion_prior_sigma
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
infection_noise_estimate <-
  function(overdispersion = FALSE,
           overdispersion_prior_mu = 0,
           overdispersion_prior_sigma = 1,
           modeldata = modeldata_init()) {
    modeldata$I_overdispersion <- overdispersion
    modeldata$I_xi_prior <- set_prior("I_xi", "normal",
      mu = overdispersion_prior_mu, sigma = overdispersion_prior_sigma
    )

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
#' @param effect_prior_mu
#' @param effect_prior_sigma
#' @param design_matrix
#'
#' @return
#' @export
#'
#' @examples
sample_effects_estimate_matrix <- function(
    design_matrix,
    effect_prior_mu = 0,
    effect_prior_sigma = 1,
    modeldata = modeldata_init()) {
  eta_prior <- set_prior(
    "eta", "normal",
    mu = effect_prior_mu, sigma = effect_prior_sigma
  )

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
#' @param effect_prior_mu
#' @param effect_prior_sigma
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
sample_effects_estimate_weekday <- function(
    effect_prior_mu = 0,
    effect_prior_sigma = 1,
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
        sample_effects_estimate_matrix(
          design_matrix, effect_prior_mu, effect_prior_sigma, modeldata
        )
    },
    required = c("meta_info$T_start_date", "meta_info$T_end_date")
  )
}

#' Title
#'
#' @param replicates
#' @param sd_prior_mu
#' @param sd_prior_sigma
#' @param replicate_sd_prior_mu
#' @param replicate_sd_prior_sigma
#' @param modeldata
#'
#' @return
#' @export
#'
#' @examples
noise_estimate <-
  function(replicates = FALSE,
           sd_prior_mu = 0,
           sd_prior_sigma = 1,
           replicate_sd_prior_mu = 0,
           replicate_sd_prior_sigma = 1,
           modeldata = modeldata_init()) {
    modeldata$replicate_noise <- replicates

    modeldata$sigma_prior <- set_prior("sigma", "truncated normal",
      mu = sd_prior_mu, sigma = sd_prior_sigma
    )

    if (modeldata$replicate_noise) {
      modeldata$tau_prior <- set_prior("tau", "truncated normal",
        mu = replicate_sd_prior_mu, sigma = replicate_sd_prior_sigma
      )
      modeldata$init$tau <- as.array(0.1)

      modeldata$init$psi <- tbe(
        rep(1e-4, modeldata$n_samples),
        "n_samples"
      )

      modeldata$checks$check_replicate_ids <- function(d) {
        if (!"replicate_ids" %in% names(d)) {
          abort(paste(
            "Replication noise can only be estimated with",
            "replicate measurements. Please specify a column `replicate_col`",
            "with replicate IDs in your data."
          ))
        }
      }
    } else {
      modeldata$tau_prior <- numeric(0)
      modeldata$init$tau <- numeric(0)
      modeldata$init$psi <- numeric(0)
    }

    modeldata$init$sigma <- 0.1 # roughly 10% coefficient of variation

    return(modeldata)
  }
