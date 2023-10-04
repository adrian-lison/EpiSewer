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
    "meta_info$length_I" = c(
      "incubation_dist_assume",
      "shedding_dist_assume",
      "residence_dist_assume"),
    "meta_info$length_R" = c(
      "incubation_dist_assume",
      "shedding_dist_assume",
      "generation_dist_assume",
      "residence_dist_assume"
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

all_components <- function() {
  components <- c(
    "concentrations",
    "droplets",
    "LOD",
    "load_per_case",
    "load_variation",
    "flows",
    "generation_dist",
    "incubation_dist",
    "shedding_dist",
    "residence_dist",
    "R",
    "seeding",
    "infection_noise",
    "sample_effects",
    "noise"
  )
  return(components)
}

#' Update meta information in modeldata based on available variables
#'
#' @description This update function is for all meta information that depends on
#'   modeldata from several functions. There are also functions which add meta
#'   information directly, in particular when this meta information cannot be
#'   solely inferred from modeldata. This function is designed such that calling
#'   it will never do harm to the modeldata object and not throw errors if
#'   something is missing in the modeldata object.
modeldata_update_metainfo <- function(modeldata) {
  if (modeldata_check(modeldata,
                      required = c("S", "D", "T"),
                      throw_error = FALSE
  )) {
    modeldata$meta_info$length_shedding <- with(modeldata, S + D + T)
  }
  if (modeldata_check(modeldata,
    required = c("L", "S", "D", "T"),
    throw_error = FALSE
  )) {
    modeldata$meta_info$length_I <- with(modeldata, L + S + D + T)
  }
  if (modeldata_check(modeldata,
    required = c("L", "S", "D", "T", "G"),
    throw_error = FALSE
  )) {
    modeldata$meta_info$length_R <- with(modeldata, L + S + D + T - G)
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

#' Construct an unspecified EpiSewer model
#'
#' @return A `modeldata` object containing data and specifications of the model
#'   to be fitted. Can be passed on to other `EpiSewer` model helpers to add
#'   further data and model specifications.
#'
#' @return The `modeldata` object also includes information about parameter
#'   initialization (`init`), meta data (`meta_info`), and checks to be
#'   performed before model fitting (`checks`).
#' @export
#'
#' @examples
#' modeldata_init()
modeldata_init <- function() {
  modeldata <- list()
  modeldata$init <- list()
  modeldata$meta_info <- list()
  modeldata$checks <- list()
  class(modeldata) <- "modeldata"
  return(modeldata)
}

#' Template function
#'
#' @description This function does not actually do anything. It only serves as a
#'   template for the documentation of other functions using inheritance of
#'   params.
#'
#' @param modeldata A `modeldata` object to which the above model specifications
#'   should be added. Default is an empty model given by [modeldata_init()]. Can
#'   also be an already partly specified model returned by other `EpiSewer`
#'   model helpers.
#'
#' @return Nothing
template_model_helpers <- function(modeldata) { }

#' Observe concentration measurements
#'
#' @description This option fits the `EpiSewer` model to pathogen concentrations
#'   measured in wastewater samples. It is suitable for different quantification
#'   methods such as qPCR or dPCR. The measured concentrations are modeled via
#'   a lognormal likelihood.
#'
#' @param data A `data.frame` with each row representing one measurement. Must
#'   have at least a column with dates and a column with concentration
#'   measurements.
#' @param composite_window Over how many days has each measured sample been
#'   collected? If 1 (default), samples represent single days. If larger than 1,
#'   samples are assumed to be equivolumetric composites over several dates. In
#'   this case, the supplied dates represent the last day included in each
#'   sample.
#' @param date_col Name of the column containing the dates.
#' @param concentration_col Name of the column containing the measured
#'   concentrations.
#' @param replicate_col Name of the column containing the replicate ID of each
#'   measurement. This is used to identify multiple measurements made of a
#'   sample from the same date. Should be `NULL` if only one measurement per
#'   date was made.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {observation types}
concentrations_observe <-
  function(data = NULL,
           composite_window = 1,
           date_col = "date",
           concentration_col = "concentration",
           replicate_col = NULL,
           modeldata = modeldata_init()) {

    if (is.null(data)) {
      data <- tryCatch(
        get_from_env("data", "measurements"),
        error = abort_f("Please supply measurement data.")
      )
    }

    if(!(composite_window%%1==0 && composite_window>0)) {
      abort("The argument `composite_window` must be a positive integer.")
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

#' Observe positive ddPCR droplet counts
#'
#' @description This option fits the `EpiSewer` model to positive droplet counts
#'   in digital droplet PCR (ddPCR) for the pathogen target of interest. This
#'   allows the use of a ddPCR-specific likelihood using a Poisson distribution.
#'   For a more generic likelihood, see [concentrations_observe()].
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @family {observation types}
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

#' Do not model a limit of detection
#'
#' @description This option drops all zero measurements from the likelihood and
#'   does not explicitly model a limit of detection.
#'
#' @details Dropping zero measurements implicitly assumes that any
#'   level of concentration can theoretically lead to a zero measurement. Since
#'   the corresponding observations are dropped, this will discard information,
#'   but not bias estimates.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {LOD models}
LOD_none <- function(modeldata = modeldata_init()) {
  modeldata$LOD <- 0
  modeldata$LOD_sharpness <- 0
  return(modeldata)
}

#' Assume a limit of detection
#'
#' @description Pathogen concentrations below a certain threshold may not be
#'   detectable and thus erroneously measured as 0. This option adjusts for a
#'   known limit of detection and includes zero measurements in the likelihood.
#'
#' @details The limit of detection is modeled using a hurdle model. In effect,
#'   zero measurements provide a signal that the concentration in the respective
#'   sample was likely below the limit of detection, but we don't know what the
#'   exact concentration was.
#'
#' @param limit Limit of detection. The concentration below which measurements
#'   will be determined as zero with substantial probability.
#' @param sharpness Sharpness of the threshold, see details.
#'
#' @details The limit of detection is highly specific to the quantification
#'   approach and protocol. It is usually established from a dedicated lab
#'   experiment.
#'
#' @details `EpiSewer` does not model a clear-cut limit of detection but rather
#'   a gradual increase in the probability of zero measurements as
#'   concentrations become smaller. This is achieved using a sigmodial curve
#'   that has its inflection point at the supplied `limit`. The `sharpness`
#'   argument determines the steepness of this curve.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso {Visualize the assumed limit of detection and its sharpness:}
#'   [plot_LOD()]
#'
#' @family {LOD models}
LOD_assume <- function(limit = NULL, sharpness = 10, modeldata = modeldata_init()) {
  if (is.null(limit)) {
    limit <- tryCatch(
      get_from_env("assumptions", "limit_of_detection"),
      error = abort_f("Please supply an assumed limit of detection.")
    )
  }
  modeldata$LOD <- limit
  modeldata$LOD_sharpness <- sharpness
  return(modeldata)
}

#' Assume the average load per case
#'
#' @description This option assumes an average total shedding load per case. In
#'   the `EpiSewer` model, this serves as a scaling factor describing how
#'   many pathogen particles are shed by the average infected individual overall
#'   and how much of this is detectable at the sampling site. This depends
#'   both on biological factors as well as on the specific sewage system.
#'
#' @param load_per_case The assumed average total shedding load per case. Must
#'   have the same unit as the numerator of the concentration unit. For example,
#'   if concentration is measured in gc/mL (gc = gene copies), then
#'   `load_per_case` should also be in gc.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso {Helper for finding a suitable load per case assumption}
#'   [suggest_load_per_case()]
load_per_case_assume <-
  function(load_per_case = NULL, modeldata = modeldata_init()) {
    if (is.null(load_per_case)) {
      load_per_case <- tryCatch(
        get_from_env("assumptions", "load_per_case"),
        error = abort_f("Please supply an assumed average shedding load per person.")
      )
    }
    modeldata$load_mean <- load_per_case
    modeldata$meta_info$load_per_case <- load_per_case
    return(modeldata)
  }

#' Do not model individual-level load variation
#'
#' @description This option assumes that there is no variation in the
#'   total load shed per case, i.e. the individual shedding load is fixed to the
#'   average shedding load.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {load variation models}
load_variation_none <- function(modeldata = modeldata_init()) {
  modeldata$load_vari <- 0
  modeldata$nu_prior <- numeric(0)
  modeldata$init$nu <- numeric(0)
  modeldata$init$zeta <- numeric(0)
  return(modeldata)
}

#' Estimate individual-level load variation
#'
#' @description This option accounts for variation in the total shedding load
#'   per case by modeling individual shedding loads as Gamma distributed with
#'   mean equal to the average `load_per_case` and a coefficient of variation to
#'   be estimated.
#'
#' @details Measurement noise and individual-level load variation might not be
#'   jointly identifiable. It may therefore be necessary to use informative
#'   priors for at least one of these parameters when using both
#'   [noise_estimate()] and [load_variation_estimate()].
#'
#' @details Also note that the accuracy of the variation model depends on
#'   the estimated number of cases to be on the right scale (which in turn
#'   depends on the assumed `load_per_case` to be roughly correct). This is
#'   because in the variation model, the population-level  coefficient of
#'   variation (CV) of shedding loads is proportional to the individual-level CV
#'   divided by the square root of the number of cases. If for example the
#'   assumed `load_per_case` is too small, the number of cases will be
#'   overestimated, which has effects in two directions:
#'   - the individual-level CV will be overestimated (especially if the prior
#'   on the individual-level CV is weak)
#'   - the population-level CV will be underestimated (especially if the prior
#'   on the individual-level CV is strong)
#'
#' @param cv_prior_mu Mean of the truncated normal prior for the coefficient of
#'   individual-level variation.
#' @param cv_prior_sigma Standard deviation of the truncated normal prior for
#'   the coefficient of individual-level variation.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {load variation models}
load_variation_estimate <- function(
    cv_prior_mu = 0.1,
    cv_prior_sigma = 0.01,
    modeldata = modeldata_init()) {
  modeldata$load_vari <- 1
  modeldata$nu_prior <- set_prior(
    "nu", "truncated normal",
    mu = cv_prior_mu, sigma = cv_prior_sigma
  )
  modeldata$init$nu <- as.array(cv_prior_mu)
  modeldata$init$zeta <- tbe(
    rep(median(modeldata$measured_concentrations, na.rm=T),
        modeldata$meta_info$length_shedding),
    c("measured_concentrations", "meta_info$length_shedding")
  )
  return(modeldata)
}

#' Assume a constant wastewater flow
#'
#' @description This option assumes a constant flow of wastewater at the
#'   sampling site. Can be used as an approximation if no regular flow
#'   measurements are available.
#'
#' @param flow_constant Daily wastewater flow volume, assumed to be the same for
#'   all days.
#'
#' @details The flow volume unit should be the same as for the concentration
#'   measurements, e.g. if concentrations are measured in gc/mL, then the flow
#'   should be in mL as well.
#'
#' @details Note that when the flow is unknown and therefore some arbitrary
#'   value (e.g. `flow_constant=1`) is assumed, this must also be accounted for
#'   in the assumed `load_per_case`. For this purpose, the function
#'   [suggest_load_per_case()] offers a `flow_constant` argument.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {flow models}
flows_assume <- function(
    flow_constant,
    modeldata = modeldata_init()) {
  modeldata <- tbc(
    "flow_data",
    {
      all_dates <-
        seq.Date(
          modeldata$meta_info$T_start_date,
          modeldata$meta_info$T_end_date,
          by = "1 day"
        )
      modeldata$flow <- rep(flow_constant, length(all_dates))
    },
    required = c("meta_info$T_start_date", "meta_info$T_end_date")
  )
  return(modeldata)
}

#' Observe wastewater flows
#'
#' @description This option accounts for daily wastewater flow volumes measured
#'   at the sampling site. The flow can change due to rainfall or industrial
#'   discharge, and directly influences pathogen concentrations in the
#'   wastewater.
#'
#' @details The flow volume unit should be the same as for the concentration
#'   measurements, e.g. if concentrations are measured in gc/mL, then the flow
#'   should be in mL as well.
#'
#' @param data A `data.frame` with each row representing one day. Must have at
#'   least a column with dates and a column with flow measurements.
#' @param date_col Name of the column containing the dates.
#' @param flow_col Name of the column containing the flows.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {flow models}
flows_observe <-
  function(data = NULL,
           date_col = "date",
           flow_col = "flow",
           modeldata = modeldata_init()) {

    if (is.null(data)) {
      data <- tryCatch(
        get_from_env("data", "flows"),
        error = abort_f("Please supply flow data.")
      )
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

#' Assume a generation time distribution
#'
#' @description This option assumes a fixed generation time distribution for the
#'   renewal model in `EpiSewer`.
#'
#' @details The generation time distribution here refers to the intrinsic
#'   distribution of the time between infection of a primary case and infection
#'   of its secondary cases. It is disease-specific and typically obtained from
#'   literature.
#'
#' @param generation_dist A numeric vector representing a discrete generation
#'   time distribution, starting with the probability for a generation time of 1
#'   day, 2 days, 3 days, and so on (a generation time of 0 days is excluded).
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso Helpers to discretize continuous probability distributions:
#'   [get_discrete_gamma()], [get_discrete_gamma_shifted()], [get_discrete_lognormal()]
generation_dist_assume <-
  function(generation_dist = NULL, modeldata = modeldata_init()) {
    if (is.null(generation_dist)) {
      generation_dist <- tryCatch(
        get_from_env("assumptions", "generation_dist"),
        error = abort_f("Please supply an assumed incubation period distribution.")
      )
    }
    modeldata$G <- length(generation_dist)
    generation_dist <- check_dist(generation_dist, "generation time distribution")
    modeldata$generation_dist <- generation_dist
    modeldata$meta_info$length_seeding <- length(generation_dist)
    return(modeldata)
  }

#' Assume an incubation period distribution
#'
#' @description This option assumes a fixed incubation period distribution for
#'   the shedding model in `EpiSewer`.
#'
#' @details `EpiSewer` uses the incubation period as a proxy for the time
#'   between infection and the start of shedding. This is because shedding load
#'   distributions in the literature are often given from symptom onset onwards.
#'   If the assumed shedding load distribution instead starts from the time of
#'   infection, then use `incubation_dist=c(1)` (i.e. no lag).
#'
#' @param incubation_dist A numeric vector representing a discrete incubation
#'   period distribution, starting with the probability for an incubation period of 0
#'   days, 1 day, 2 days, and so on.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso Helpers to discretize continuous probability distributions:
#'   [get_discrete_gamma()], [get_discrete_lognormal()]
incubation_dist_assume <-
  function(incubation_dist = NULL, modeldata = modeldata_init()) {
    if (is.null(incubation_dist)) {
      incubation_dist <- tryCatch(
        get_from_env("assumptions", "incubation_dist"),
        error = abort_f("Please supply an assumed incubation period distribution.")
      )
    }
    modeldata$L <- length(incubation_dist) - 1
    incubation_dist <- check_dist(incubation_dist, "incubation period distribution")
    modeldata$incubation_dist <- incubation_dist
    return(modeldata)
  }

#' Assume a shedding load distribution
#'
#' @description This option assumes a fixed shedding load distribution for the
#'   shedding model in `EpiSewer`.
#'
#' @param shedding_dist A numeric vector representing a discrete shedding load
#'   distribution, with elements describing the share of load shed 0 days, 1
#'   day, 2 days, and so on after the start of shedding.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso Helpers to discretize continuous probability distributions:
#'   [get_discrete_gamma()], [get_discrete_lognormal()]
shedding_dist_assume <-
  function(shedding_dist = NULL, modeldata = modeldata_init()) {
    if (is.null(shedding_dist)) {
      shedding_dist <- tryCatch(
        get_from_env("assumptions", "shedding_dist"),
        error = abort_f("Please supply an assumed shedding load distribution.")
      )
    }
    modeldata$S <- length(shedding_dist) - 1
    shedding_dist <- check_dist(shedding_dist, "shedding load distribution")
    modeldata$shedding_dist <- shedding_dist
    return(modeldata)
  }

#' Assume a sewer residence time distribution
#'
#' @description This option assumes a fixed residence time distribution for
#'   pathogen particles. By default, `EpiSewer` assumes that particles arrive at
#'   the sampling site within the day of shedding. However, for larger sewage
#'   systems, particles may travel longer than a day depending on where and
#'   when they were shed into the wastewater.
#'
#' @param residence_dist A numeric vector representing a discrete residence time
#'   distribution, with elements describing the share of load that takes 0 days,
#'   1 day, 2 days, and so on to arrive at the sampling site.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @examples
#' # Particles arrive within the same day
#' residence_dist_assume(residence_dist = c(1))
#'
#' # Particles always arrive after one day
#' residence_dist_assume(residence_dist = c(0, 1))
#'
#' 1/4 of particles only arrives after one day
#' residence_dist_assume(residence_dist = c(0.75, 0.25))
residence_dist_assume <-
  function(residence_dist = NULL, modeldata = modeldata_init()) {
    if (is.null(residence_dist)) {
      residence_dist <- tryCatch(
        get_from_env("assumptions", "residence_dist"),
        error = abort_f("Please supply an assumed residence time distribution.")
      )
    }
    modeldata$D <- length(residence_dist) - 1
    residence_dist <- check_dist(residence_dist, "residence time distribution")
    modeldata$residence_dist <- residence_dist
    return(modeldata)
  }

#' Smooth Rt estimates via exponential smoothing
#'
#' @description This option estimates the effective reproduction number over
#'   time using exponential smoothing. It implements Holt's linear trend method
#'   with dampening through an innovations state space model with a level,
#'   trend, and dampening component.
#'
#' @param level_prior_mu Prior (mean) on the initial level of Rt.
#' @param level_prior_sigma Prior (standard deviation) on the initial level of
#'   Rt.
#' @param trend_prior_mu Prior (mean) on the initial trend of Rt.
#' @param trend_prior_sigma Prior (standard deviation) on the initial trend of
#'   Rt.
#' @param sd_prior_mu Prior (mean) on the standard deviation of the innovations.
#' @param sd_prior_sigma Prior (standard deviation) on the standard deviation of
#'   the innovations.
#' @param smooth_fixed Fixed value for the smoothing parameter
#'   (alpha). If NULL (default), the smoothing parameter is estimated, otherwise
#'   it will be fixed to the given value.
#' @param smooth_prior_shapes Prior (Beta shape 1 and 2) on the smoothing
#'   parameter.
#' @param trend_smooth_fixed Fixed value for the trend smoothing parameter
#'   (beta). If NULL (default), the trend smoothing parameter is estimated,
#'   otherwise it will be fixed to the given value.
#' @param trend_smooth_prior_shapes Prior (Beta shape 1 and 2) on the trend
#'   smoothing parameter.
#' @param dampen_fixed Fixed value for the dampening parameter (phi). If NULL
#'   (default), the dampening parameter is estimated, otherwise it will be fixed
#'   to the given value.
#' @param dampen_prior_shapes Prior (Beta shape 1 and 2) on the dampening
#'   parameter.
#' @param differenced If `FALSE` (default), exponential smoothing is applied to
#'   the absolute Rt time series. If `TRUE`, it is instead applied to the
#'   differenced time series. This makes the level become the trend, and the
#'   trend become the curvature.
#' @param noncentered If `TRUE` (default), a non-centered parameterization is
#'   used to model the innovations in the state space process (for better
#'   sampling efficiency).
#'
#' @details The innovations state space model consists of three components: a
#'   level, a trend, and a dampening component.
#' - The level is smoothed based on the levels from earlier time steps,
#'   with exponentially decaying weights, as controlled by a smoothing parameter
#'   (often called alpha). Note that *smaller* values of `alpha` indicate
#'   *stronger* smoothing. In particular, `alpha = 1` means that only the last
#'   level is used.
#' - The trend is smoothed based on the trends from earlier time steps,
#'   with exponentially decaying weights, as controlled by a trend smoothing
#'   parameter (often called beta). Note that *smaller* values of `beta`
#'   indicate *stronger* smoothing. In particular, `beta = 1` means that only
#'   the last trend is used.
#' - The dampening determines how long a previous trend continues into the
#'   future before it levels of to a stationary time series. The strength of
#'   dampening is controlled by a dampening parameter (often called phi).
#'   Note that *smaller* values of `phi` indicate *stronger* dampening.
#'   In particular, `phi = 1` means no dampening. Values below `phi = 0.8`
#'   are seldom in practice as the dampening becomes very strong.
#'
#' @details Often, `alpha`, `beta`, and `phi` are jointly unidentifiable. It may
#'   therefore be necessary to fix at least one of the parameters (typically
#'   `phi`) or supply strong priors.
#'
#' @details Note that the smoothness of retrospective Rt estimates is often more
#'   influenced by the prior on the standard deviation of innovations than the
#'   smoothing and trend smoothing parameters. The smoothing parameters mostly
#'   have an influence on estimates close to the present / date of estimation,
#'   when limited data signal is available. Here, the standard deviation of the
#'   innovations influences how uncertain Rt estimates are close to the present.
#'
#' @details The priors of this component have the following functional form:
#' - initial level of Rt: `Normal`
#' - initial trend of Rt: `Normal`
#' - standard deviation of innovations: `Truncated normal`
#' - smoothing parameter: `Beta`
#' - trend smoothing parameter: `Beta`
#' - dampening parameter: `Beta`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
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

#' Smooth Rt estimates via a random walk
#'
#' @description This option estimates the effective reproduction number over
#'   time using a random walk.
#'
#' @param intercept_prior_mu Prior (mean) on the intercept of the random walk.
#' @param intercept_prior_sigma Prior (standard deviation) on the intercept of
#'   the random walk.
#' @param sd_prior_mu Prior (mean) on the standard deviation of the random walk.
#' @param sd_prior_sigma Prior (standard deviation) on the standard deviation of
#'   the random walk.
#' @param differenced If `FALSE` (default), the random walk is applied to the
#'   absolute Rt time series. If `TRUE`, it is instead applied to the
#'   differenced time series, i.e. now the trend is modeled as a random walk.
#' @param noncentered If `TRUE` (default), a non-centered parameterization is
#'   used to model the innovations of the random walk (for better sampling
#'   efficiency).
#'
#' @details The smoothness of Rt estimates is influenced by the prior on the
#'   standard deviation of the random walk. It also influences the uncertainty
#'   of Rt estimates towards the present / date of estimation, when limited
#'   data signal is available. The prior on the intercept of the random walk
#'   should reflect your expectation of Rt at the beginning of the time series.
#'   If estimating from the start of an epidemic, you might want to use a prior
#'   with mean > 1 for the intercept.
#'
#' @details The priors of this component have the following functional form:
#' - intercept of the random walk: `Normal`
#' - standard deviation of the random walk: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
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

#' Smooth Rt estimates via smoothing splines
#'
#' @description This option estimates the effective reproduction number using
#'   regularized smoothing splines on a B-spline basis.
#'
#' @param knot_distance Distance between spline breakpoints (knots) in days
#'   (default is 1, i.e. a knot on each day). Placing knots further apart
#'   increases the smoothness of Rt estimates and can speed up model fitting. If
#'   knots are too far apart however, the estimates can become inaccurate.
#' @param spline_degree Degree of the spline polynomials (default is 3 for cubic
#'   splines).
#' @param coef_intercept_prior_mu Prior (mean) on the intercept of the log
#'   random walk over spline coefficients.
#' @param coef_intercept_prior_sigma Prior (standard deviation) on the intercept
#'   of the log random walk over spline coefficients.
#' @param coef_sd_prior_mu Prior (mean) on the standard deviation of the random
#'   walk over spline coefficients.
#' @param coef_sd_prior_sigma Prior (standard deviation) on the standard
#'   deviation of the random walk over spline coefficients.
#'
#' @details `EpiSewer` uses a random walk on the B-spline coefficients for
#'   regularization. This allows to use small knot distances without obtaining
#'   extremely wiggly Rt estimates. The random walk is modeled on the
#'   logarithmic scale (i.e. a geometric random walk):
#'   - A prior on the log intercept with mean 0 roughly equals a prior on
#'   the unit intercept with mean 1. The prior on the intercept should reflect
#'   your expectation of Rt at the beginning of the time series. If estimating
#'   from the start of an epidemic, you might want to use a prior with
#'   log(mean) > 0 (i.e. mean > 1) for the intercept.
#'   - The prior on the standard deviation should be interpreted in terms of
#'   exponential multiplicative changes. For example, a standard deviation of
#'   0.2 on the log scale roughly allows for multiplicative changes between
#'   exp(-0.4) and exp(0.4) on the unit scale (using the 2 sigma rule).
#'
#' @details The smoothness of the Rt estimates is influenced both by the knot
#'   distance and by the standard deviation of the random walk on coefficients.
#'   The latter also influences the uncertainty of Rt estimates towards the
#'   present / date of estimation, when limited data signal is available.
#'
#' @details The priors of this component have the following functional form:
#' - intercept of the log random walk over spline coefficients: `Normal`
#' - standard deviation of the log random walk over spline coefficients: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
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

#' Estimate seeding infections
#'
#' @description This option estimates initial infections at the start of the
#'   modeled time period, when the renewal model cannot be applied yet. It uses
#'   a random walk to model these seeding infections.
#'
#' @param intercept_prior_mu Prior (mean) on the intercept of the seeding random
#'   walk. If NULL (default), this is set to a crude empirical estimate of the
#'   number of cases (see details).
#' @param intercept_prior_sigma Prior (standard deviation) on the intercept of
#'   the seeding random walk.
#' @param sd_prior_mu Prior (mean) on the standard deviation of the seeding
#'   random walk.
#' @param sd_prior_sigma Prior (standard deviation) on the standard deviation of
#'   the seeding random walk.
#'
#' @details The seeding phase has the length of the maximum generation time
#'   (during this time, the renewal model cannot be applied). Traditionally,
#'   seeding refers to the first few (potentially imported) infections of an
#'   epidemic, but depending on what time period the model is fitted to, this
#'   may also cover a different phase with higher incidence levels.
#'
#' @details If `intercept_prior_mu` is not specified by the user, `EpiSewer`
#'   will set it to a rough estimate of the number of cases using the supplied
#'   wastewater measurements and shedding assumptions. We note that this is a
#'   violation of Bayesian principles (data must not be used to inform priors) -
#'   but a neglectable one, since it only ensures that the seeding is modeled on
#'   the right order of magnitude and does not have relevant impacts on later Rt
#'   estimates.
#'
#' @details The priors of this component have the following functional form:
#' - intercept of the seeding random walk: `Normal`
#' - standard deviation of the seeding random walk: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
seeding_estimate <- function(
    intercept_prior_mu = NULL,
    intercept_prior_sigma = 1,
    sd_prior_mu = 0.05,
    sd_prior_sigma = 0.025,
    modeldata = modeldata_init()) {
  if (!is.null(intercept_prior_mu)) {
    modeldata$iota_log_ar_start_prior <- set_prior(
      "iota_log_ar_start",
      "normal",
      mu = intercept_prior_mu,
      sigma = intercept_prior_sigma
    )
  } else {
    modeldata$iota_log_ar_start_prior <- tbe(
      set_prior(
        "iota_log_ar_start",
        "normal (mu based on crude empirical estimate of cases)",
        mu = log(modeldata$meta_info$initial_cases_crude),
        sigma = intercept_prior_sigma
      ),
      "meta_info$initial_cases_crude"
    )
  }

  modeldata$iota_log_ar_sd_prior <- set_prior(
    "iota_log_ar_sd",
    "truncated normal",
    mu = sd_prior_mu,
    sigma = sd_prior_sigma
  )

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

#' Do not model infection noise
#'
#' This option does not model noise in the infection process and instead
#' implements a deterministic renewal model.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {infection noise models}
infection_noise_none <- function(modeldata = modeldata_init()) {
  modeldata$I_sample <- FALSE
  modeldata$I_overdispersion <- FALSE
  modeldata$init$I <- numeric(0)
  modeldata$init$I_log <- numeric(0)
  return(modeldata)
}

#' Estimate infection noise
#'
#' @description This option estimates noise in the infection process, i.e.
#'   implements a stochastic renewal model. This allows for variation in the
#'   number of new infections generated at each time step, which can often speed
#'   up model fitting.
#'
#' @param overdispersion If `FALSE` (default) new infections are modeled as
#'   Poisson distributed. If `TRUE`, new infections are modeled as Negative
#'   Binomial distributed.
#' @param overdispersion_prior_mu Prior (mean) on the overdispersion parameter
#'   of the Negative Binomial.
#' @param overdispersion_prior_sigma Prior (standard deviation) on the
#'   overdispersion parameter of the Negative Binomial.
#'
#' @details The priors of this component have the following functional form:
#' - overdispersion parameter of the Negative Binomial: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {infection noise models}
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

#' Do not model sample effects
#'
#' @description This option does not model effects of sample covariates on the
#'   concentrations.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {sample effect models}
sample_effects_none <- function(modeldata = modeldata_init()) {
  modeldata$K <- 0
  modeldata$X <- numeric(0)
  modeldata$eta_prior <- numeric(0)
  modeldata$init$eta <- numeric(0)
  return(modeldata)
}

#' Estimate sample effects using a design matrix
#'
#' @description This option uses a linear regression model to estimate effects
#'   of sample covariates on the concentration. Concentrations can be influenced
#'   by sampling-related external factors, for example the time between sampling
#'   and shipping to the lab (age-of-sample effect), or different sampling or
#'   storage methods.
#'
#' @param effect_prior_mu Prior (mean) on the regression coefficients.
#' @param effect_prior_sigma Prior (standard deviation) on the regression
#'   coefficients.
#' @param design_matrix A design matrix with different covariates potentially
#'   influencing sample concentration. The matrix must have one row for each
#'   modeled day. See [stats::model.matrix()] for construction of design
#'   matrices.
#'
#' @details `EpiSewer` will fit a fixed-effects linear model, random effects are
#'   currently not supported.
#'
#' @details The priors of this component have the following functional form:
#' - regression coefficients: `Normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {sample effect models}
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
            "days (from earliest to latest date, accounting for the composite window length),",
            "but design matrix for sample date effects has",
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

#' Estimate weekday sample effects
#'
#' @description This option uses a linear regression model to estimate sample
#'   weekday effects on the concentration. Concentrations can be influenced by
#'   the time between sampling and shipping to the lab (age-of-sample effect),
#'   and if shipment follows a weekly batch scheme, the sampling weekday is a
#'   good proxy for the age at shipment.
#'
#' @param effect_prior_mu Prior (mean) on the regression coefficients.
#' @param effect_prior_sigma Prior (standard deviation) on the regression
#'   coefficients.
#'
#' @details Effects are estimated for weekdays Monday - Saturday, with Sunday as
#'   the baseline. `EpiSewer` will fit a fixed-effects linear model, random
#'   effects are currently not supported.
#'
#' @details The priors of this component have the following functional form:
#' - regression coefficients: `Normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {sample effect models}
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

#' Estimate measurement noise
#'
#' @description This option estimates the unexplained variation in wastewater
#'   measurements. If multiple measurements (replicates) per sample are
#'   provided, `EpiSewer` can also explicitly model variation before the
#'   replication stage.
#'
#' @param replicates Should replicates be used to explicitly model variation
#'   before the replication stage?
#' @param cv_prior_mu Prior (mean) on the coefficient of variation of
#'   concentration measurements.
#' @param cv_prior_sigma Prior (standard deviation) on the coefficient of
#'   variation of concentration measurements.
#' @param pre_replicate_cv_prior_mu Prior (mean) on the coefficient of variation
#'   of concentrations before the replication stage.
#' @param pre_replicate_cv_prior_sigma Prior (standard deviation) on the
#'   coefficient of variation of concentrations before the replication stage.
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements: `Truncated normal`
#' - coefficient of variation of concentration before the replication stage:
#' `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
noise_estimate <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           modeldata = modeldata_init()) {
    modeldata$pr_noise <- replicates

    modeldata$cv_prior <- set_prior("cv", "truncated normal",
      mu = cv_prior_mu, sigma = cv_prior_sigma
    )

    if (modeldata$pr_noise) {
      modeldata$tau_prior <- set_prior("tau", "truncated normal",
        mu = pre_replicate_cv_prior_mu, sigma = pre_replicate_cv_prior_sigma
      )
      modeldata$init$tau <- as.array(pre_replicate_cv_prior_mu)

      modeldata$init$psi <- tbe(
        rep(1e-4, modeldata$n_samples),
        "n_samples"
      )

      modeldata$checks$check_replicate_ids <- function(d) {
        if (!"replicate_ids" %in% names(d)) {
          abort(paste(
            "Variation before the replication stage can only be estimated with",
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

    modeldata$init$cv <- 0.1 # 10% coefficient of variation

    return(modeldata)
  }
