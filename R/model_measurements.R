#' Model the measurement process
#'
#' @description This module function is used to specify the components of the
#'  `measurements` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param concentrations Concentration measurements from wastewater samples.
#' Modeling options:
#' `r component_functions_("concentrations")`
#' @param noise Measurement noise due to unexplained variation in sampling
#' and lab analysis. Modeling options:
#' `r component_functions_("noise")`
#' @param LOD Limit of detection. Concentrations below a certain threshold may
#' not be detectable and thus erroneously measured as 0. `EpiSewer` can adjust
#' for the limit of detection using a hurdle model. Modeling options:
#' `r component_functions_("LOD")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `measurements` module.
#' @export
model_measurements <- function(
    concentrations = concentrations_observe(),
    noise = noise_estimate(),
    LOD = LOD_none()) {
  verify_is_modeldata(concentrations, "concentrations")
  verify_is_modeldata(noise, "noise")
  verify_is_modeldata(LOD, "LOD")
  return(modeldata_combine(concentrations, noise, LOD))
}

#' Observe concentration measurements
#'
#' @description This option fits the `EpiSewer` model to pathogen concentrations
#'   measured in wastewater samples. It is suitable for different quantification
#'   methods such as qPCR or dPCR. The measured concentrations are modeled via a
#'   lognormal likelihood.
#'
#' @param measurements A `data.frame` with each row representing one
#'   measurement. Must have at least a column with dates and a column with
#'   concentration measurements.
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
  function(measurements = NULL,
           composite_window = 1,
           date_col = "date",
           concentration_col = "concentration",
           replicate_col = NULL,
           modeldata = modeldata_init()) {
    if (!(composite_window %% 1 == 0 && composite_window > 0)) {
      rlang::abort(
        "The argument `composite_window` must be a positive integer."
      )
    }

    modeldata <- tbp("concentrations_observe",
      {
        required_data_cols <- c(date_col, concentration_col, replicate_col)
        if (!all(required_data_cols %in% names(measurements))) {
          rlang::abort(
            c(paste(
              "The following columns must be present",
              "in the provided measurements `data.frame`:",
              paste(required_data_cols, collapse = ", ")
            ),
            paste("Please adjust the `data.frame` or specify the right column",
                  "names via the `_col` arguments of this function."))
          )
        }

        if (nrow(measurements)==0) {
          rlang::abort("The provided measurements `data.frame` is empty.")
        }

        measurements[[date_col]] <- as.Date(measurements[[date_col]])

        if (is.null(replicate_col)) {
          if (any(duplicated(measurements[[date_col]]))) {
            rlang::abort(
              paste(
                "Duplicated dates found in measurements `data.frame`.",
                "If your data contains replicate measurements, please add a",
                "column with replicate IDs to your `data.frame` and provide its",
                "name via the `replicate_col` argument."
              )
            )
          }
        }

        modeldata$T <-
          as.integer(max(measurements[[date_col]]) -
            min(measurements[[date_col]]) + composite_window)

        modeldata$.metainfo$T_start_date <-
          min(measurements[[date_col]]) - composite_window + 1
        modeldata$.metainfo$T_end_date <- max(measurements[[date_col]])

        modeldata$w <- composite_window
        modeldata$.metainfo$composite_window <- composite_window


        measured <- !is.na(measurements[[concentration_col]])
        modeldata$n_measured <- sum(measured)
        modeldata$n_samples <- length(
          unique(measurements[[date_col]][measured])
          )
        measured_dates <- as.integer(
          measurements[[date_col]][measured] -
            modeldata$.metainfo$T_start_date + 1
        )
        modeldata$sample_to_date <- sort(unique(measured_dates))
        modeldata$measure_to_sample <- sapply(
          measured_dates,
          function(x) which(x == modeldata$sample_to_date)[[1]]
        )
        modeldata$measured_concentrations <- as.numeric(
          measurements[[concentration_col]][measured]
        )
        modeldata$.metainfo$measured_dates <- as.Date(
          measurements[[date_col]][measured]
        )

        if (!is.null(replicate_col)) {
          modeldata$replicate_ids <- as.integer(
            measurements[[replicate_col]][measured]
          )
        }

        return(modeldata)
      },
      required_data = "measurements",
      modeldata = modeldata
    )

    if (composite_window != 1) {
      .str_details <- c(composite_window = composite_window)
    } else {
      .str_details <- c()
    }
    modeldata$.str$measurements[["concentrations"]] <- list(
      concentrations_observe = .str_details
    )

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
    rlang::abort(paste(
      "Specification of measurements via ddPCR droplet count",
      "is not implemented yet."
    ))
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

    modeldata$nu_upsilon_prior <- set_prior("nu_upsilon", "truncated normal",
      mu = cv_prior_mu, sigma = cv_prior_sigma
    )

    if (modeldata$pr_noise) {
      modeldata$nu_psi_prior <- set_prior(
        "nu_psi", "truncated normal",
        mu = pre_replicate_cv_prior_mu, sigma = pre_replicate_cv_prior_sigma
      )
      modeldata$.init$nu_psi <- as.array(pre_replicate_cv_prior_mu)

      modeldata$.init$psi <- tbe(
        rep(1e-4, modeldata$n_samples),
        "n_samples"
      )

      modeldata$.checks$check_replicate_ids <- function(d) {
        if (!"replicate_ids" %in% names(d)) {
          rlang::abort(paste(
            "Variation before the replication stage can only be estimated with",
            "replicate measurements. Please specify a column `replicate_col`",
            "with replicate IDs in your data."
          ))
        }
      }
    } else {
      modeldata$nu_psi_prior <- numeric(0)
      modeldata$.init$nu_psi <- numeric(0)
      modeldata$.init$psi <- numeric(0)
    }

    modeldata$.init$nu_upsilon <- 0.1 # 10% coefficient of variation

    if (replicates == TRUE) {
      .str_details <- c(replicates = replicates)
    } else {
      .str_details <- c()
    }
    modeldata$.str$measurements[["noise"]] <- list(
      noise_estimate = .str_details
    )

    return(modeldata)
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

  modeldata$.str$measurements[["LOD"]] <- list(
    LOD_none = c()
  )

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
LOD_assume <- function(limit = NULL, sharpness = 10,
                       modeldata = modeldata_init()) {

  if (sharpness<=0) {
    rlang::abort("Sharpness parameter must be greater than zero.")
  }
  limit_of_detection <- limit
  modeldata <- tbp("LOD_assume",
    {
      modeldata$LOD <- limit
      modeldata$LOD_sharpness <- sharpness
      return(modeldata)
    },
    required_assumptions = "limit_of_detection",
    modeldata = modeldata
  )

  modeldata$.str$measurements[["LOD"]] <- list(
    LOD_assume = c()
  )

  return(modeldata)
}
