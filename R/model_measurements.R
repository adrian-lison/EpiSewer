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
#'   Currently supported are "normal" (default and recommended), and
#'   "lognormal". Distributions are truncated below zero if necessary.
#' @param date_col Name of the column containing the dates.
#' @param concentration_col Name of the column containing the measured
#'   concentrations.
#' @param replicate_col Name of the column containing the replicate ID of each
#'   measurement. This is used to identify multiple measurements made of a
#'   sample from the same date. Should be `NULL` if only one measurement per
#'   date was made.
#' @param n_averaged The number of replicates over which the measurements have
#'   been averaged. This is typically used as an alternative to providing
#'   several replicates per sample (i.e. the concentration provided in the
#'   `measurements` `data.frame` is the arithmetic mean of several replicates).
#'   Can be either a single number (it is then assumed that the number of
#'   averaged replicates is the same for each observation) or a vector (one
#'   value for each observation).
#' @param n_averaged_col Name of the column in the `measurements` data.frame
#'   containing the number of replicates over which the measurements have been
#'   averaged. This is an alternative to specifying `n_averaged`.
#' @param total_droplets_col Name of the column in the `measurements`
#'   data.frame containing the total number of droplets in the ddPCR reaction
#'   of each measurement. Only applies to concentration measurements obtain via
#'   ddPCR. Can be used by the [noise_estimate_ddPCR()] and
#'   [LOD_estimate_ddPCR()] modeling components. Note that this is really the
#'   *total* number of droplets, not just the number of positive droplets.
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
           n_averaged = 1,
           n_averaged_col = NULL,
           total_droplets_col = NULL,
           modeldata = modeldata_init()) {
    if (!(composite_window %% 1 == 0 && composite_window > 0)) {
      cli::cli_abort(
        "The argument `composite_window` must be a positive integer."
      )
    }

    if (!distribution %in% c("normal", "lognormal")) {
      cli::cli_abort(
        c(
          "Only the following distributions are supported:",
          "normal",
          "lognormal"
        )
      )
    }

    modeldata <- tbp("concentrations_observe",
      {
        required_data_cols <- c(
          date_col, concentration_col, replicate_col,
          n_averaged_col, total_droplets_col
          )
        if (!all(required_data_cols %in% names(measurements))) {
          cli::cli_abort(
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
          cli::cli_abort("The provided measurements `data.frame` is empty.")
        }

        measurements[[date_col]] <- as.Date(measurements[[date_col]])

        if (is.null(replicate_col)) {
          if (any(duplicated(measurements[[date_col]]))) {
            cli::cli_abort(
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

        if (!is.null(n_averaged_col)){
          if (!n_averaged_col %in% names(measurements)) {
            cli::cli_abort(paste0(
                "The column `", n_averaged_col, "` does not exist in the ",
                "measurements `data.frame`."
              ))
          }
          modeldata$n_averaged <- as.numeric(measurements[[n_averaged_col]][measured])
          if (any(is.na(modeldata$n_averaged))) {
            cli::cli_abort(paste0(
                "The column `", n_averaged_col, "` contains missing ",
                "values for some of the observed measurements."
              ))
          }
        } else if (length(n_averaged) == 1) {
          modeldata$n_averaged <- rep(n_averaged, modeldata$n_measured)
        } else if (length(n_averaged) == modeldata$n_measured) {
          modeldata$n_averaged <- n_averaged
        } else {
          cli::cli_abort(paste(
              "The length of `n_averaged` must be either 1 or equal to the",
              "number of samples."
            ))
        }

        if (!is.null(total_droplets_col)) {
          modeldata$ddPCR_total_droplets <- as.integer(
            measurements[[total_droplets_col]][measured]
          )
        } else {
          modeldata$ddPCR_total_droplets <- rep(0, modeldata$n_measured)
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
    cli::cli_abort(paste(
      "Specification of measurements via ddPCR droplet count",
      "is not implemented yet."
    ))
  }

#' Estimate measurement noise (internal helper function)
#'
#' @description This option estimates the unexplained variation in wastewater
#'   measurements. If multiple measurements (replicates) per sample are
#'   provided, `EpiSewer` can also explicitly model variation before the
#'   replication stage.
#'
#' @description This helper function is called from [noise_estimate()] and
#'   [noise_estimate_ddPCR()]. [noise_estimate()] is a constant coefficient of
#'   variation model, [noise_estimate_ddPCR()] is a noise model specialized for
#'   digital droplet PCR (`cv_type = "ddPCR"`), which may however also work with
#'   other quantification methods such as qPCR.
#'
#' @param replicates Should replicates be used to explicitly model variation
#'   before the replication stage?
#' @param cv_prior_mu Prior (mean) on the coefficient of variation of
#'   concentration measurements. Note that when `replicates=TRUE`, this is only
#'   the CV after the replication stage (see details for more explanation).
#' @param cv_prior_sigma Prior (standard deviation) on the coefficient of
#'   variation of concentration measurements.
#' @param cv_type One out of "constant" (default), "constant_var", or "ddPCR".
#'   If "constant", the coefficient of variation is estimated as a
#'   constant/single parameter for all observations. If "ddPCR", the coefficient
#'   of variation is modeled as a function of the expected concentration
#'   according to the statistical properties of ddPCR. In particular, this model
#'   predicts a higher coefficient of variation at smaller concentrations, which
#'   often leads to a better model fit. If "constant_var", not the coefficient
#'   of variation but the variance of measurements is modeled as constant. This
#'   is usually a misspecification and is only supported for comparison
#'   purposes.
#' @param ddPCR_prior_droplets_mu Prior (mean) on the total number of droplets
#'   in the ddPCR reaction.
#' @param ddPCR_prior_droplets_sigma Prior (standard deviation) on the total
#'   number of droplets in the ddPCR reaction. If this is set to zero, the total
#'   number of droplets will be fixed to the prior mean and not estimated.
#' @param ddPCR_droplets_observe If TRUE, the number of total droplets is taken
#'   from the supplied measurements `data.frame`. This requires that the
#'   argument `total_droplets_col` is specified in [concentrations_observe()].
#' @param ddPCR_droplet_variation_prior_mu Prior (mean) on the coefficient of
#'   variation of the total number of droplets in the ddPCR reaction. Usually,
#'   the maximum number of partitions possible for a given dPCR chip is not
#'   reached, i.e. a certain number of partitions is lost. This loss varies
#'   between PCR runs, and is modeled as log-normal distributed in EpiSewer.
#' @param ddPCR_droplet_variation_prior_sigma Prior (standard deviation) on the
#'   coefficient of variation of the total number of droplets in the ddPCR
#'   reaction.
#' @param ddPCR_prior_scaling_mu Prior (mean) on the concentration scaling
#'   factor for the ddPCR reaction. The concentration scaling factor is the
#'   droplet volume, scaled by the dilution of the wastewater in the ddPCR
#'   reaction. See details for further explanation.
#' @param ddPCR_prior_scaling_sigma Prior (standard deviation) on the
#'   concentration scaling factor for the ddPCR reaction. If this is set to
#'   zero, the concentration scaling factor will be fixed to the prior mean and
#'   not estimated.
#' @param pre_replicate_cv_prior_mu Prior (mean) on the coefficient of variation
#'   of concentrations *before* the replication stage.
#' @param pre_replicate_cv_prior_sigma Prior (standard deviation) on the
#'   coefficient of variation of concentrations *before* the replication stage.
#'
#' @param prePCR_noise_type The parametric distribution to assume for noise
#'   before the PCR assay. Currently supported are "log-normal" and "gamma". The
#'   choice of the parametric distribution typically makes no relevant
#'   difference for the noise model, but can make a relevant difference for the
#'   LOD model if [LOD_estimate_ddPCR()] is used.
#'
#' @param use_taylor_approx If TRUE (default), a Taylor expansion approximation
#'   is used to estimate the CV of measurements under pre-PCR noise. The
#'   approximation is very accurate, unless concentrations are extremely high
#'   (so high that the quality of the measurements from ddPCR would anyway be
#'   questionable).
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
noise_estimate_ <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           cv_type = "constant",
           ddPCR_prior_droplets_mu = NULL,
           ddPCR_prior_droplets_sigma = NULL,
           ddPCR_droplets_observe = NULL,
           ddPCR_droplet_variation_prior_mu = NULL,
           ddPCR_droplet_variation_prior_sigma = NULL,
           ddPCR_prior_scaling_mu = NULL,
           ddPCR_prior_scaling_sigma = NULL,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           prePCR_noise_type = "log-normal",
           use_taylor_approx = TRUE,
           modeldata = modeldata_init()) {
    modeldata$pr_noise <- replicates

    modeldata$nu_upsilon_a_prior <- set_prior(
      "nu_upsilon_a", "truncated normal",
      mu = cv_prior_mu, sigma = cv_prior_sigma
    )
    modeldata$.init$nu_upsilon_a <- 0.1 # 10% coefficient of variation

    if (cv_type == "constant") {
      modeldata$ddPCR_droplets_observe <- FALSE
      modeldata$cv_type <- 0
      modeldata$nu_upsilon_b_mu_prior <- numeric(0)
      modeldata$nu_upsilon_b_cv_prior <- numeric(0)
      modeldata$.init$nu_upsilon_b_mu <- numeric(0)
      modeldata$.init$nu_upsilon_b_cv <- numeric(0)
      modeldata$.init$nu_upsilon_b_noise_raw <- numeric(0)
      modeldata$nu_upsilon_c_prior <- numeric(0)
      modeldata$nu_upsilon_c_fixed <- -1 # dummy value
      modeldata$.init$nu_upsilon_c <- numeric(0)
      modeldata$cv_pre_type <- numeric(0)
      modeldata$cv_pre_approx_taylor <- numeric(0)
    } else if (cv_type == "ddPCR") {
      modeldata$cv_type <- 1

      if (ddPCR_droplets_observe) {
        modeldata$ddPCR_droplets_observe <- TRUE
        modeldata$nu_upsilon_b_mu_prior <- set_prior(
          "nu_upsilon_b_mu", "dummy prior", mu = 0, sigma = 0
        )
        modeldata$nu_upsilon_b_cv_prior <- numeric(0)
        modeldata$.init$nu_upsilon_b_mu <- numeric(0)
        modeldata$.init$nu_upsilon_b_cv <- numeric(0)
        modeldata$.init$nu_upsilon_b_noise_raw <- numeric(0)
        modeldata$.checks$check_total_droplets_col <- function(d) {
          if (!"ddPCR_total_droplets" %in% names(d)) {
            cli::cli_abort(paste0(
              "You specified `ddPCR_droplets_observe = TRUE`, which requires ",
              "a column with the number of total droplets in the PCR for each ",
              "sample in your data. Please specify such a column via the",
              "`total_droplets_col` argument in ",
              cli_help("concentrations_observe"), "."
            ))
          }
        }
      } else {
        modeldata$ddPCR_droplets_observe <- FALSE

        modeldata$nu_upsilon_b_mu_prior <- set_prior(
          "nu_upsilon_b_mu", "truncated normal",
          mu = ddPCR_prior_droplets_mu * 1e-4, # scaling by 1e-4 for numerical reasons
          sigma = ddPCR_prior_droplets_sigma * 1e-4
        )
        modeldata$.init$nu_upsilon_b_mu <- init_from_normal_prior(
          modeldata$nu_upsilon_b_mu_prior
          )

        modeldata$nu_upsilon_b_cv_prior <- set_prior(
          "nu_upsilon_b_cv", "truncated normal",
          mu = ddPCR_droplet_variation_prior_mu,
          sigma = ddPCR_droplet_variation_prior_sigma
        )
        modeldata$.init$nu_upsilon_b_cv <- as.array(0.01)
        modeldata$.init$nu_upsilon_b_noise_raw <- tbe(
          rep(0, modeldata$n_measured), "n_measured"
        )
      }

      modeldata$nu_upsilon_c_prior <- set_prior(
        "nu_upsilon_c", "truncated normal",
        mu = ddPCR_prior_scaling_mu * 1e+5, # scaling by 1e+5 for numerical reasons
        sigma = ddPCR_prior_scaling_sigma * 1e+5
      )
      modeldata$.init$nu_upsilon_c <- init_from_normal_prior(
        modeldata$nu_upsilon_c_prior
        )

      if (prePCR_noise_type == "gamma") {
        modeldata$cv_pre_type <- 0
      } else if (prePCR_noise_type == "log-normal") {
        modeldata$cv_pre_type <- 1
      } else {
        cli::cli_abort(paste0(
            "`prePCR_noise_type = ", prePCR_noise_type, "` not supported.",
            "Available options: 'gamma', `log-normal`."
          ))
      }
      modeldata$cv_pre_approx_taylor <- use_taylor_approx

    } else if (cv_type == "constant_var") {
      modeldata$ddPCR_droplets_observe <- FALSE
      modeldata$cv_type <- 2
      modeldata$nu_upsilon_b_mu_prior <- numeric(0)
      modeldata$nu_upsilon_b_cv_prior <- numeric(0)
      modeldata$.init$nu_upsilon_b_mu <- numeric(0)
      modeldata$.init$nu_upsilon_b_cv <- numeric(0)
      modeldata$.init$nu_upsilon_b_noise_raw <- numeric(0)
      modeldata$nu_upsilon_c_prior <- numeric(0)
      modeldata$nu_upsilon_c_fixed <- -1 # dummy value
      modeldata$.init$nu_upsilon_c <- numeric(0)
      modeldata$cv_pre_type <- numeric(0)
      modeldata$cv_pre_approx_taylor <- numeric(0)
    } else {
      cli::cli_abort(
        paste0(
          "Noise type `", cv_type, "` not supported. Available options: ",
          "'constant', `ddPCR`, `constant_var`."
        ),
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
          cli::cli_abort(paste(
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

#' Estimate measurement noise
#'
#' @description This option estimates the unexplained variation in wastewater
#'   measurements using a constant coefficient of variation model.
#'
#' @description If multiple measurements (replicates) per sample are provided,
#'   `EpiSewer` can also explicitly model variation before the replication
#'   stage.
#'
#' @description For a non-constant coefficient of variation model, see
#'   [noise_estimate_ddPCR()].
#'
#' @details When `replicates=TRUE`, two coefficients of variation are estimated:
#' - the CV before the replication stage (see `pre_replicate_cv_prior_mu`)
#' - the CV after the replication stage (see `cv_prior_mu`)
#'
#' The meaning of these CV estimates depends on the type of replicates. If the
#' replicates are biological replicates (i.e. independently processed), then
#' `cv` estimates the noise in the preprocessing and PCR, and
#' `pre_replicate_cv` estimates the noise from anything before preprocessing
#' (e.g. sampling noise and all other unexplained variation). In contrast, if
#' the replicates are technical replicates (i.e. several PCR runs of the same
#' preprocessed sample), then `cv` estimates only the PCR noise, and
#' `pre_replicate_cv` estimates all other noise (including preprocessing
#' noise.)
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements (`cv`): `Truncated normal`
#' - coefficient of variation of concentration before the replication stage (`pre_replicate_cv`):
#'   `Truncated normal`
#'
#' @inheritParams noise_estimate_
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family [noise models]
noise_estimate <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           modeldata = modeldata_init()) {
    return(noise_estimate_(
      replicates = replicates,
      cv_prior_mu = cv_prior_mu,
      cv_prior_sigma = cv_prior_sigma,
      cv_type = "constant",
      pre_replicate_cv_prior_mu = pre_replicate_cv_prior_mu,
      pre_replicate_cv_prior_sigma = pre_replicate_cv_prior_sigma,
      modeldata = modeldata
      ))
  }

#' Estimate measurement noise for digital droplet PCR data
#'
#' @description This option estimates the unexplained variation in wastewater
#'   measurements using a coefficient of variation model specialised for digital
#'   droplet PCR. Specifically, the coefficient of variation is modeled as a
#'   function of the expected concentration according to the statistical
#'   properties of ddPCR.
#'
#' @description This model predicts a higher coefficient of variation at smaller
#'   concentrations, which often leads to a better model fit, even for
#'   measurements from other quantification methods such as qPCR.
#'
#' @description If multiple measurements (replicates) per sample are provided,
#'   `EpiSewer` can also explicitly model variation before the replication
#'   stage.
#'
#' @param cv_prior_mu Prior (mean) on the coefficient of variation of
#'   concentration measurements. Note that in contrast to using
#'   [noise_estimate()], this does *not* include the technical noise of the
#'   digital PCR. This is because the dPCR noise is explicitly modeled (using
#'   the number of partitions and conversion factor). Moreover, when
#'   `replicates=TRUE`, this is only the CV after the replication stage (see
#'   details for more explanation).
#'
#' @details The concentration scaling factor (see `ddPCR_prior_scaling_mu`,
#'   `ddPCR_prior_scaling_sigma`) is the droplet volume,
#'   scaled by the dilution of the wastewater in the ddPCR reaction. The
#'   dilution accounts for all extraction and preparation steps. For example, if
#'   the droplet volume is 4.5e-7 mL and the dilution of the wastewater is 100:3
#'   (i.e. 100 gc/mL in the original wastewater sample correspond to 3 gc/mL in
#'   the PCR reaction), then the overall scaling factor is 4.5e-7 * 100 / 3 =
#'   1.5e-5.
#'
#' @details When `replicates=TRUE`, two coefficients of variation are estimated:
#' - the CV before the replication stage (see `pre_replicate_cv_prior_mu`)
#' - the CV after the replication stage (see `cv_prior_mu`)
#'
#' The meaning of these CV estimates depends on the type of replicates. If the
#' replicates are biological replicates (i.e. independently processed), then
#' `cv` estimates the noise in the preprocessing before the PCR, and
#' `pre_replicate_cv` estimates the noise from anything before preprocessing
#' (e.g. sampling noise and all other unexplained variation). In contrast, if
#' the replicates are technical replicates (i.e. several PCR runs of the same
#' preprocessed sample), then `cv` estimates only unexplained PCR noise
#' (should be close to zero), and `pre_replicate_cv` estimates all other noise
#' (including preprocessing noise.)
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements (`cv`): `Truncated normal`
#' - mean number of droplets in ddPCR: `Truncated normal`
#' - coefficient of variation of number of droplets in ddPCR: `Truncated normal`
#' - concentration scaling factor of ddPCR: `Truncated normal`
#' - coefficient of variation of concentration before the replication stage (`pre_replicate_cv`):
#'   `Truncated normal`
#'
#' @inheritParams noise_estimate_
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family [noise models]
#' @seealso [LOD_estimate_ddPCR] for a limit of detection model specialised for
#'   ddPCR.
noise_estimate_ddPCR <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           ddPCR_prior_droplets_mu = 20000,
           ddPCR_prior_droplets_sigma = 5000,
           ddPCR_droplets_observe = FALSE,
           ddPCR_droplet_variation_prior_mu = 0,
           ddPCR_droplet_variation_prior_sigma = 0.05,
           ddPCR_prior_scaling_mu = 1.5e-5,
           ddPCR_prior_scaling_sigma = 0.5e-5,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           prePCR_noise_type = "log-normal",
           use_taylor_approx = TRUE,
           modeldata = modeldata_init()) {
    return(noise_estimate_(
      replicates = replicates,
      cv_prior_mu = cv_prior_mu,
      cv_prior_sigma = cv_prior_sigma,
      cv_type = "ddPCR",
      ddPCR_prior_droplets_mu = ddPCR_prior_droplets_mu,
      ddPCR_prior_droplets_sigma = ddPCR_prior_droplets_sigma,
      ddPCR_droplets_observe = ddPCR_droplets_observe,
      ddPCR_droplet_variation_prior_mu = ddPCR_droplet_variation_prior_mu,
      ddPCR_droplet_variation_prior_sigma = ddPCR_droplet_variation_prior_sigma,
      ddPCR_prior_scaling_mu = ddPCR_prior_scaling_mu,
      ddPCR_prior_scaling_sigma = ddPCR_prior_scaling_sigma,
      pre_replicate_cv_prior_mu = pre_replicate_cv_prior_mu,
      pre_replicate_cv_prior_sigma = pre_replicate_cv_prior_sigma,
      prePCR_noise_type = prePCR_noise_type,
      use_taylor_approx = use_taylor_approx,
      modeldata = modeldata
    ))
  }

#' Estimate measurement noise with constant variance
#'
#' @description This option estimates the unexplained variation in wastewater
#'   measurements using a constant variance model. This is usually a
#'   misspecification and is only supported for comparison purposes.
#'
#' @description For a constant coefficient of variation model, see
#'   [noise_estimate()], and for a non-constant coefficient of variation model,
#'   see [noise_estimate_ddPCR()].
#'
#' @description If multiple measurements (replicates) per sample are provided,
#'   `EpiSewer` can also explicitly model variation before the replication
#'   stage.
#'
#' @details Note that although this model keeps the variance constant, the prior
#'   for the measurement noise is still in terms of the (average) coefficient of
#'   variation (CV). This makes prior specification easier since it the CV is
#'   unitless.
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements: `Truncated normal`
#' - coefficient of variation of concentration before the replication stage:
#'   `Truncated normal`
#'
#' @inheritParams noise_estimate_
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family [noise models]
noise_estimate_constant_var <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           warn = TRUE,
           modeldata = modeldata_init()) {
    if (warn) {
      cli::cli_warn(paste(
      "You have specified",
      cli_help("noise_estimate_constant_var"),
      "as the model component for measurement noise.",
      "Note that modeling a constant variance",
      "is likely a model misspecification and should only be used for ",
      "comparison purposes with better models like",
      cli_help("noise_estimate"), "or", cli_help("noise_estimate_ddPCR"), ".",
      "You can specify",
      "{.code noise_estimate_constant_var(warn=TRUE)} to disable this warning."
      ))
    }
    return(noise_estimate_(
      replicates = replicates,
      cv_prior_mu = cv_prior_mu,
      cv_prior_sigma = cv_prior_sigma,
      cv_type = "constant_var",
      pre_replicate_cv_prior_mu = pre_replicate_cv_prior_mu,
      pre_replicate_cv_prior_sigma = pre_replicate_cv_prior_sigma,
      modeldata = modeldata
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
  modeldata$LOD_model <- 0
  modeldata$LOD_scale <- numeric(0)
  modeldata$LOD_drop_prob <- 0

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
#' @param drop_prob Probability for non-detection below which likelihood
#'   contributions of observed concentrations are dropped from LOD model. This
#'   avoids numerical issues of the LOD model at high concentrations (very small
#'   non-detection probabilities) that can otherwise affect sampling speed.
#'   Since these likelihood contributions will be virtually zero for almost all
#'   samples anyway, parameter estimates are practically not affected.
#'
#' @details The limit of detection is specific to the quantification approach
#'   and protocol. It is usually established from a dedicated lab experiment
#'   (serial dilution experiment). It his here assumed that this experiment did
#'   not cover a large fraction of the preprocessing noise to find an optimal
#'   configuration for the exponential model.
#'
#' @details If used together with [noise_estimate_ddPCR()], EpiSewer will also
#'   model the effect of pre-PCR noise on the LOD. This means that the modeled
#'   LOD could be slightly higher than specified under `limit`, depending on the
#'   estimated pre-PCR noise.
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
                       drop_prob = 1e-10,
                       modeldata = modeldata_init()) {

  if (!LOD_type %in% c("exponential", "ddPCR")) { # "ddPCR" is synonym
    cli::cli_abort(
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

  modeldata$LOD_drop_prob <- drop_prob

  modeldata$.str$measurements[["LOD"]] <- list(
    LOD_assume = c()
  )

  return(modeldata)
}

#' Estimate a limit of detection model for digital droplet PCR data
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
#' @param drop_prob Probability for non-detection below which likelihood
#'   contributions of observed concentrations are dropped from LOD model. This
#'   avoids numerical issues of the LOD model at high concentrations (very small
#'   non-detection probabilities) that can otherwise affect sampling speed.
#'   Since these likelihood contributions will be virtually zero for almost all
#'   samples anyway, parameter estimates are practically not affected.
#'
#' @details The limit of detection is modeled using a hurdle model. The model
#'   uses the number of droplets in the ddPCR reaction and the concentration
#'   scaling factor as defined and estimated by [noise_estimate_ddPCR()]. It can
#'   therefore only be used together with `noise = noise_estimate_ddPCR()` in
#'   [model_measurements()].
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @family {LOD models}
LOD_estimate_ddPCR <- function(drop_prob = 1e-10, modeldata = modeldata_init()) {

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
      'To use ',
      cli_help("LOD_estimate_ddPCR"), ', you must specify noise = ',
      cli_help("noise_estimate_ddPCR"), ' in ', cli_help("model_measurements"),
      "."
      ),
    modeldata = modeldata
    )

  modeldata$LOD_drop_prob <- drop_prob

  modeldata$.str$measurements[["LOD"]] <- list(
    LOD_estimate_ddPCR = c()
  )

  return(modeldata)
}
