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
#' @param distribution Parametric distribution for concentration measurements.
#'   Currently supported are normal (default and recommended), and lognormal.
#'   Distributions are truncated below zero if necessary.
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
           distribution = "normal",
           date_col = "date",
           concentration_col = "concentration",
           replicate_col = NULL,
           modeldata = modeldata_init()) {
    if (!(composite_window %% 1 == 0 && composite_window > 0)) {
      rlang::abort(
        "The argument `composite_window` must be a positive integer."
      )
    }

    if (!distribution %in% c("normal", "lognormal")) {
      rlang::abort(
        c(
          "Only the following distributions are supported:",
          "normal",
          "lognormal"
        )
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

    modeldata$obs_dist = switch(
      distribution,
      "normal" = 1,
      "lognormal" = 2,
      -1
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
#' @description Aside from a constant coefficient of variation model, a
#'   noise model specialized for digital droplet PCR (`cv_type = "ddPCR"`) is
#'   available, which may however also work with other quantification
#'   methods such as qPCR.
#'
#' @param replicates Should replicates be used to explicitly model variation
#'   before the replication stage?
#' @param cv_prior_mu Prior (mean) on the coefficient of variation of
#'   concentration measurements.
#' @param cv_prior_sigma Prior (standard deviation) on the coefficient of
#'   variation of concentration measurements.
#' @param cv_type One out of "constant" (default), or "ddPCR". If ddPCR, the
#'   coefficient of variation is modeled as a function of the expected
#'   concentration according to the statistical properties of ddPCR. In
#'   particular, this model predicts a higher coefficient of variation at
#'   smaller concentrations, which often leads to a better model fit.
#' @param ddPCR_prior_droplets_mu Prior (mean) on the number of droplets in the
#'   ddPCR reaction.
#' @param ddPCR_prior_droplets_sigma Prior (standard deviation) on the number of
#'   droplets in the ddPCR reaction.
#' @param ddPCR_droplets_fixed If TRUE (default), the number of droplets is
#'   fixed to the prior mean and not estimated. This is recommended if no
#'   replicates are available.
#' @param ddPCR_prior_scaling_mu Prior (mean) on the concentration scaling
#'   factor for the ddPCR reaction. The concentration scaling factor is the
#'   droplet volume, scaled by the dilution of the wastewater in the ddPCR
#'   reaction. See details for further explanation.
#' @param ddPCR_prior_scaling_sigma Prior (standard deviation) on the
#'   concentration scaling factor for the ddPCR reaction.
#' @param ddPCR_scaling_fixed If TRUE, the concentration scaling factor is fixed
#'   to the prior mean and not estimated.
#' @param pre_replicate_cv_prior_mu Prior (mean) on the coefficient of variation
#'   of concentrations before the replication stage.
#' @param pre_replicate_cv_prior_sigma Prior (standard deviation) on the
#'   coefficient of variation of concentrations before the replication stage.
#'
#' @details The concentration scaling factor (see `ddPCR_prior_scaling_mu`,
#'   `ddPCR_prior_scaling_sigma`, `ddPCR_scaling_fixed`) is the droplet volume,
#'   scaled by the dilution of the wastewater in the ddPCR reaction. The
#'   dilution accounts for all extraction and preparation steps. For example, if
#'   the droplet volume is 4.5e-4 and the dilution of the wastewater is 30 (i.e.
#'   30 gc/mL in the original wastewater sample correspond to 1 gc/mL in the PCR
#'   reaction), then the overall scaling factor is 4.5e-4 / 30 = 1.5e-5.
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements: `Truncated normal`
#' - coefficient of variation of concentration before the replication stage:
#'   `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
noise_estimate <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           cv_type = "constant",
           ddPCR_prior_droplets_mu = 20000,
           ddPCR_prior_droplets_sigma = 5000,
           ddPCR_droplets_fixed = TRUE,
           ddPCR_prior_scaling_mu = 1.5e-5,
           ddPCR_prior_scaling_sigma = 0.5e-5,
           ddPCR_scaling_fixed = FALSE,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           modeldata = modeldata_init()) {
    modeldata$pr_noise <- replicates

    modeldata$nu_upsilon_a_prior <- set_prior("nu_upsilon_a", "truncated normal",
      mu = cv_prior_mu, sigma = cv_prior_sigma
    )
    modeldata$.init$nu_upsilon_a <- 0.1 # 10% coefficient of variation

    if (cv_type == "constant") {
      modeldata$cv_type <- 0
      modeldata$nu_upsilon_b_prior <- numeric(0)
      modeldata$nu_upsilon_b_fixed <- -1 # dummy value
      modeldata$.init$nu_upsilon_b <- numeric(0)
      modeldata$nu_upsilon_c_prior <- numeric(0)
      modeldata$nu_upsilon_c_fixed <- -1 # dummy value
      modeldata$.init$nu_upsilon_c <- numeric(0)
    } else if (cv_type == "ddPCR") {
      modeldata$cv_type <- 1
      if (ddPCR_droplets_fixed) {
        modeldata$nu_upsilon_b_fixed <- ddPCR_prior_droplets_mu * 1e-4
        modeldata$nu_upsilon_b_prior <- numeric(0)
        modeldata$.init$nu_upsilon_b <- numeric(0)
      } else {
        modeldata$nu_upsilon_b_fixed <- -1 # dummy value
        modeldata$nu_upsilon_b_prior <- set_prior(
          "nu_upsilon_b", "truncated normal",
          mu = ddPCR_prior_droplets_mu * 1e-4, # scaling by 1e-4 for numerical reasons
          sigma = ddPCR_prior_droplets_sigma * 1e-4
        )
        modeldata$.init$nu_upsilon_b <- as.array(1)
      }

      if (ddPCR_scaling_fixed) {
        modeldata$nu_upsilon_c_fixed <- ddPCR_prior_scaling_mu * 1e+5
        modeldata$nu_upsilon_c_prior <- numeric(0)
        modeldata$.init$nu_upsilon_c <- numeric(0)
      } else {
        modeldata$nu_upsilon_c_fixed <- -1 # dummy value
        modeldata$nu_upsilon_c_prior <- set_prior(
          "nu_upsilon_c", "truncated normal",
          mu = ddPCR_prior_scaling_mu * 1e+5, # scaling by 1e+5 for numerical reasons
          sigma = ddPCR_prior_scaling_sigma * 1e+5
        )
        modeldata$.init$nu_upsilon_c <- as.array(1)
      }

    } else {
      rlang::abort(
      paste0("Noise type `", cv_type, "` not supported. Available options: 'constant', `ddPCR`."),
      )
    }

    if (modeldata$pr_noise) {
      modeldata$nu_psi_prior <- set_prior(
        "nu_psi", "truncated normal",
        mu = pre_replicate_cv_prior_mu, sigma = pre_replicate_cv_prior_sigma
      )
      modeldata$.init$nu_psi <- as.array(0.1)

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
  modeldata$LOD_model <- 0
  modeldata$LOD_scale <- numeric(0)

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
#' @param limit Limit of detection. The concentration below which the pathogen
#'   cannot be detected with sufficient probability, i.e. the measurement may be
#'   zero although the pathogen is present in the sample.
#' @param prob What desired probability of detection does the limit refer to?
#'   Default is 95% (0.95): This means that the provided `limit` is the smallest
#'   concentration at which the pathogen can still be detected with over 95%
#'   probability.
#' @param LOD_type The type of LOD model used. Currently, only "exponential" is
#'   supported. This models an exponentially decreasing probability of zero
#'   measurements / non-detection as a function of concentration. The
#'   exponential model can be derived from the statistical properties of ddPCR,
#'   but should also work well for other quantification methods such as qPCR.
#'
#' @details The limit of detection is specific to the quantification approach
#'   and protocol. It is usually established from a dedicated lab experiment.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso {Visualize the assumed LOD as a function of concentration:}
#'   [plot_LOD()]
#'
#' @family {LOD models}
LOD_assume <- function(limit = NULL, prob = 0.95, LOD_type = "exponential",
                       modeldata = modeldata_init()) {

  if (!LOD_type %in% c("exponential", "ddPCR")) { # "ddPCR" is synonym
    rlang::abort(
      'Currently, only LOD_type = "exponential" is supported.'
    )
  }

  limit_of_detection <- limit
  modeldata <- tbp("LOD_assume",
    {
      modeldata$LOD_model <- 1
      modeldata$LOD_scale <- -log(1-prob)/limit
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

#' Estimate a limit of detection model for ddPCR data
#'
#' @description Pathogen concentrations below a certain threshold may not be
#'   detectable and thus erroneously measured as 0. This option adjusts for a
#'   limit of detection based on the statistical properties of digital droplet
#'   PCR (ddPCR) and includes zero measurements in the likelihood.
#'
#' @description In effect, zero measurements provide a signal that the
#'   concentration in the respective sample was likely below the limit of
#'   detection, but we don't know what the exact concentration was.
#'
#' @details The limit of detection is modeled using a hurdle model. The model
#'   uses the number of droplets in the ddPCR reaction and the concentration
#'   scaling factor as defined and estimated by [noise_estimate()]. It is
#'   therefore necessary to specify `cv_type = "ddPCR"` in [noise_estimate()].
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @family {LOD models}
LOD_estimate_ddPCR <- function(modeldata = modeldata_init()) {

  modeldata <- tbc("LOD_estimate_ddPCR",
    {
      modeldata$LOD_model <- 2
      modeldata$LOD_scale <- numeric(0)
      return(modeldata)
    },
    required = c(
      "cv_type"
      ),
    required_values = c(1),
    advice = paste0(
      'To use {.run [LOD_estimate_ddPCR()](EpiSewer::LOD_estimate_ddPCR()}, ',
      'you must specify cv_type = "ddPCR" in ',
      '{.run [noise_estimate()](EpiSewer::noise_estimate()}.'
      ),
    modeldata = modeldata
    )

  modeldata$.str$measurements[["LOD"]] <- list(
    LOD_estimate_ddPCR = c()
  )

  return(modeldata)
}
