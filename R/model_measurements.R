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
#' @family {module functions}
model_measurements <- function(
    concentrations = concentrations_observe(),
    noise = noise_estimate(),
    LOD = LOD_none()) {
  verify_is_modeldata(concentrations, "concentrations")
  verify_is_modeldata(noise, "noise")
  verify_is_modeldata(LOD, "LOD")
  return(modeldata_combine(concentrations, noise, LOD))
}

#' Observe measurements
#'
#' @description This helper function is called from other measurement modeling
#'   functions. Currently supported are [concentrations_observe()] and
#'   [concentrations_observe_partitions()].
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
#'   Currently supported are "gamma" (default and recommended), "log-normal",
#'   "truncated normal", and "normal". The "truncated normal" and "normal"
#'   options are not recommended for use in practice.
#' @param date_col Name of the column containing the dates.
#' @param concentration_col Name of the column containing the measured
#'   concentrations.
#' @param positive_partitions_col Name of the column in the `measurements`
#'   data.frame containing the number of positive partitions (e.g. positive
#'   droplets for ddPCR) in the dPCR reaction of each measurement. If several
#'   technical replicates are used, this should be the AVERAGE number of
#'   positive partitions per replicate.
#' @param replicate_col Name of the column containing the replicate ID of each
#'   measurement. This is used to identify multiple measurements made of a
#'   sample from the same date. Should be `NULL` if only one measurement per
#'   date was made.
#' @param n_averaged The number of technical replicates (i.e. repeated PCR runs)
#'   used for each sample or biological replicate. The concentration provided in
#'   the `measurements` `data.frame` is assumed to be the average or pooled
#'   estimate from several technical replicates. Can be either a single number
#'   (it is then assumed that the number of averaged replicates is the same for
#'   each observation) or a vector (one value for each observation).
#' @param n_averaged_col Name of the column in the `measurements` data.frame
#'   containing the number of technical replicates over which the measurements
#'   have been averaged/pooled. This is an alternative to specifying
#'   `n_averaged`.
#' @param total_partitions_col Name of the column in the `measurements`
#'   data.frame containing the number of total partitions (e.g. droplets for
#'   ddPCR) in the dPCR reaction of each measurement. If several technical
#'   replicates are used, this should be the AVERAGE number of valid partitions
#'   per replicate. Only applies when modeling concentration measurements via
#'   the dPCR-specific noise model. Can be used by the [noise_estimate_dPCR()]
#'   and [LOD_estimate_dPCR()] modeling components. Note that this is really the
#'   number of *valid* partitions, not the number of positive partitions.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
measurements_observe_ <- function(
    measurements = NULL,
    observation_type = c("concentrations", "partitions"),
    composite_window = 1,
    distribution = "gamma",
    date_col = "date",
    concentration_col = NULL,
    positive_partitions_col = NULL,
    replicate_col = NULL,
    n_averaged = 1,
    n_averaged_col = NULL,
    total_partitions_col = NULL,
    modeldata = modeldata_init()) {

    observation_type <- rlang::arg_match(observation_type)

    if (!(composite_window %% 1 == 0 && composite_window > 0)) {
      cli::cli_abort(
        "The argument `composite_window` must be a positive integer."
      )
    }

    modeldata <- tbp("measurements_observe",
      {
        if (observation_type == "concentrations") {
          required_data_cols <- list(
            date_col, concentration_col, replicate_col,
            n_averaged_col, total_partitions_col
          )
          data_col_names <- c(
            "date", "concentration", "replicate_id",
            "n_averaged", "total_partitions"
          )[!sapply(required_data_cols, is.null)]
        } else if (observation_type == "partitions") {
          required_data_cols <- list(
            date_col, positive_partitions_col, concentration_col,
            replicate_col, n_averaged_col, total_partitions_col
          )
          data_col_names <- c(
            "date", "positive_partitions", "concentration",
            "replicate_id", "n_averaged", "total_partitions"
          )[!sapply(required_data_cols, is.null)]
        }

        required_data_cols <- purrr::list_c(required_data_cols)
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
        measurements = as.data.table(measurements)[, .SD, .SDcols = required_data_cols]
        measurements <- setnames(measurements, old = required_data_cols, new = data_col_names)

        if (observation_type == "concentrations") {
          modeldata$.metainfo$measurements_cols <- list(
            date_col = date_col,
            concentration_col = concentration_col,
            replicate_col = replicate_col,
            n_averaged_col = n_averaged_col,
            total_partitions_col = total_partitions_col
          )
        } else if (observation_type == "partitions") {
          modeldata$.metainfo$measurements_cols <- list(
            date_col = date_col,
            positive_partitions_col = positive_partitions_col,
            concentration_col = concentration_col,
            replicate_col = replicate_col,
            n_averaged_col = n_averaged_col,
            total_partitions_col = total_partitions_col
          )
        }

        # housekeeping
        measurements <- measurements[!is.na(concentration) & !is.na(date), ]
        measurements[, date := as.Date(date)]
        measurements[, concentration := as.numeric(concentration)]
        if (observation_type == "partitions") {
          measurements[, positive_partitions := as.numeric(positive_partitions)]
        }

        if (nrow(measurements)==0) {
          cli::cli_abort(
            "The provided measurements `data.frame` contains no measurements."
            )
        }

        if (is.null(replicate_col)) {
          if (any(duplicated(measurements[["date"]]))) {
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

        modeldata$T <- as.integer(
          measurements[, max(date) - min(date) + composite_window]
          )

        modeldata$.metainfo$T_start_date <- measurements[
          , min(date) - composite_window + 1
          ]
        modeldata$.metainfo$T_end_date <- measurements[, max(date)]

        modeldata$w <- composite_window
        modeldata$.metainfo$composite_window <- composite_window

        modeldata$n_measured <- nrow(measurements)
        modeldata$n_samples <- length(unique(measurements[["date"]]))
        measured_dates <- as.integer(
          measurements[["date"]] - modeldata$.metainfo$T_start_date + 1
        )
        modeldata$sample_to_date <- sort(unique(measured_dates))
        modeldata$measure_to_sample <- sapply(measured_dates, function(x) {
          which(x == modeldata$sample_to_date)[[1]]
        })
        modeldata$.metainfo$measured_dates <- measurements[["date"]]

        modeldata$measured_concentrations <- measurements[["concentration"]]
        if (observation_type == "concentrations") {
          modeldata$positive_partitions <- numeric(0)
          modeldata$.init$concentration_with_noise_raw <- numeric(0)
        } else if (observation_type == "partitions") {
          modeldata$positive_partitions <- measurements[["positive_partitions"]]
          modeldata$.init$concentration_with_noise_raw <- rep(1, modeldata$n_measured)
        }

        if (!is.null(replicate_col)) {
          modeldata$replicate_ids <- as.integer(measurements[["replicate_id"]])
        }

        if (!is.null(n_averaged_col)) {
          modeldata$n_averaged <- as.numeric(measurements[["n_averaged"]])
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

        # total valid partitions in PCR run
        modeldata <- tbc("dPCR_total_partitions", {
          if (modeldata$total_partitions_observe) {
            if (!is.null(total_partitions_col)) {
              modeldata$dPCR_total_partitions <- as.integer(
                measurements[["total_partitions"]]
              )
            } else {
              cli::cli_abort(paste0(
                "You specified `total_partitions_observe = TRUE`, but this ",
                "requires a column with the total number of partitions in the ",
                "dPCR run. Please specify a column with the total number of ",
                "partitions via the `total_partitions_col` ",
                "argument in ", cli_help("concentrations_observe"), " or ",
                cli_help("concentrations_observe_partitions"), "."
              ))
            }
          } else {
            if (observation_type == "partitions") {
              cli::cli_abort(paste0(
                "You specified `total_partitions_observe = FALSE`, but ",
                "total partition counts are required for the ",
                "`concentrations_observe_partitions` model component. ",
                "Please specify `total_partitions_observe = TRUE` ",
                "in ", cli_help("noise_estimate_dPCR"),
                " and provide ",
                "a column with the total number of partitions via the ",
                "`total_partitions_col` argument in ",
                cli_help("concentrations_observe_partitions"), "."
              ))
            }
            if (!is.null(total_partitions_col)) {
              cli::cli_inform(c("i" = paste0(
                "Note: Your data contains a column with the number of total ",
                "partitions in the dPCR."), "!" = paste0("However, you specified ",
                 "total_partitions_observe = FALSE, so this column is currently ignored."
                )))
            }
            modeldata$dPCR_total_partitions <- numeric(0)
          }
        }, required = c("total_partitions_observe"), modeldata = modeldata)

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

    if (observation_type == "concentrations") {
      modeldata$.str$measurements[["concentrations"]] <- list(
        concentrations_observe = .str_details
      )
    } else if (observation_type == "partitions") {
      modeldata$.str$measurements[["concentrations"]] <- list(
        concentrations_observe_partitions = .str_details
      )
    }

    available_distributions <- c(
      "gamma" = 0,
      "log-normal" = 1,
      "truncated normal" = 2,
      "normal" = 3,
      "binomial" = 4
    )

    distribution_aliases <- c(
      "lnorm" = "log-normal",
      "lognormal" = "log-normal",
      "truncated_normal" = "truncated normal",
      "truncated-normal" = "truncated normal",
      "norm" = "normal"
    )

    distribution <- stringr::str_to_lower(distribution)
    distribution <- do.call(
      switch, c(distribution, as.list(c(distribution_aliases,distribution)))
    )

    if (!distribution %in% names(available_distributions)) {
      cli::cli_abort(
        c(
          "Only the following distributions are supported:",
          setNames(names(available_distributions), rep("*", length(available_distributions)))
        )
      )
    }

    if (distribution == "truncated normal") {
      cli::cli_inform(c(
        "!" = paste(
          "The Truncated Normal distribution is supported for model",
          "comparison purposes -",
          "but it is not recommended in practice due to misspecification, i.e.",
          "biased mean at low concentrations.",
          "Consider using a Gamma or Log-Normal distribution instead."
      )))
    }
    if (distribution == "normal") {
      cli::cli_inform(c(
        "!" = paste(
          "The Normal distribution is supported for model",
          "comparison purposes -",
          "but it is not recommended in practice due to misspecification, i.e.",
          "prediction of negative concentrations.",
          "Consider using a Gamma or Log-Normal distribution instead."
      )))
    }
    if (distribution == "binomial" && observation_type != "partitions") {
      cli::cli_abort(paste0(
        "The Binomial distribution is only supported for positive partitions ",
        "counts in dPCR assays. Please use `concentrations_observe_partitions()` ",
        "to specify this distribution type."
      ))
    }

    modeldata$obs_dist = do.call(
      switch, c(distribution, as.list(c(available_distributions,-1)))
      )

    return(modeldata)
}

#' Observe concentration measurements
#'
#' @description This option fits the `EpiSewer` model to pathogen concentrations
#'   measured in wastewater samples. It is suitable for different quantification
#'   methods such as qPCR or dPCR. By default, the measured concentrations are
#'   modeled via a gamma likelihood.
#'
#' @inheritParams template_model_helpers
#' @inheritParams measurements_observe_
#' @inherit modeldata_init return
#'
#' @export
#' @family {observation types}
concentrations_observe <- function(
    measurements = NULL,
    composite_window = 1,
    distribution = "gamma",
    date_col = "date",
    concentration_col,
    replicate_col = NULL,
    n_averaged = 1,
    n_averaged_col = NULL,
    total_partitions_col = NULL,
    modeldata = modeldata_init()) {
  return(measurements_observe_(
    measurements = measurements,
    observation_type = "concentrations",
    composite_window = composite_window,
    distribution = distribution,
    date_col = date_col,
    concentration_col = concentration_col,
    positive_partitions_col = NULL,
    replicate_col = replicate_col,
    n_averaged = n_averaged,
    n_averaged_col = n_averaged_col,
    total_partitions_col = total_partitions_col,
    modeldata = modeldata))
}

#' Observe positive dPCR partition counts
#'
#' @description This option fits the `EpiSewer` model to positive partition
#'   counts in digital PCR (dPCR), e.g. positive droplets in ddPCR. This allows
#'   the use of a dPCR-specific likelihood using a Binomial model for the number
#'   of positive partitions observed. For a more generic likelihood, see
#'   [concentrations_observe()].
#'
#' @inheritParams template_model_helpers
#' @inheritParams measurements_observe_
#' @inherit modeldata_init return
#'
#' @export
#' @family {observation types}
concentrations_observe_partitions <- function(
    measurements = NULL,
    composite_window = 1,
    date_col = "date",
    concentration_col,
    positive_partitions_col,
    replicate_col = NULL,
    n_averaged = 1,
    n_averaged_col = NULL,
    total_partitions_col,
    modeldata = modeldata_init()) {
      return(measurements_observe_(
        measurements = measurements,
        observation_type = "partitions",
        composite_window = composite_window,
        distribution = "binomial",
        date_col = date_col,
        concentration_col = concentration_col,
        positive_partitions_col = positive_partitions_col,
        replicate_col = replicate_col,
        n_averaged = n_averaged,
        n_averaged_col = n_averaged_col,
        total_partitions_col = total_partitions_col,
        modeldata = modeldata))
  }

#' Estimate measurement noise (internal helper function)
#'
#' @description This option estimates the unexplained variation in wastewater
#'   measurements. If multiple measurements (replicates) per sample are
#'   provided, `EpiSewer` can also explicitly model variation before the
#'   replication stage.
#'
#' @description This helper function is called from other noise modeling
#'   functions. [noise_estimate()] is a constant coefficient of variation model,
#'   [noise_estimate_dPCR()] is a noise model specialized for digital PCR
#'   (`cv_type = "dPCR"`), which may however also work with other quantification
#'   methods such as qPCR, and [noise_estimate_constant_var()] is a constant
#'   variance model.
#'
#' @param replicates Should replicates be used to explicitly model variation
#'   before the replication stage?
#' @param cv_prior_mu Prior (mean) on the coefficient of variation of
#'   concentration measurements. Note that when `replicates=TRUE`, this is only
#'   the CV after the replication stage (see details for more explanation).
#' @param cv_prior_sigma Prior (standard deviation) on the coefficient of
#'   variation of concentration measurements.
#' @param cv_type One out of "constant" (default), "constant_var", or "dPCR". If
#'   "constant", the coefficient of variation is estimated as a constant/single
#'   parameter for all observations. If "dPCR", the coefficient of variation is
#'   modeled as a function of the expected concentration according to the
#'   statistical properties of dPCR. In particular, this model predicts a higher
#'   coefficient of variation at smaller concentrations, which often leads to a
#'   better model fit. If "constant_var", not the coefficient of variation but
#'   the variance of measurements is modeled as constant. This is usually a
#'   misspecification and is only supported for comparison purposes.
#' @param total_partitions_observe If TRUE, the total number of partitions is
#'   taken from the supplied measurements `data.frame`. This requires that the
#'   argument `total_partitions_col` is specified in [concentrations_observe()].
#' @param max_partitions_prior_lower Prior (5% quantile) for the maximum total
#'   number of dPCR partitions. This is usually defined by the manufacturer of
#'   the dPCR system/chip used, which supports a certain maximum number of
#'   partitions. If you know the exact dPCR system and its maximum partition
#'   number, you can set both `max_partitions_prior_lower` and
#'   `max_partitions_prior_upper` to this value. Otherwise, this prior can be
#'   used to set a broad lower and upper bound for the maximum number of
#'   partitions, to reflect a range of popular dPCR systems/chips.
#' @param max_partitions_prior_upper Prior (95% quantile) for the maximum total
#'   number of dPCR partitions (see `max_partitions_prior_lower` for details.)
#' @param partition_loss_mean_prior_lower Prior (5% quantile) for the mean
#'   relative partition loss. A certain proportion of partitions in a dPCR run
#'   is typically invalid and discarded from the concentration estimate. This
#'   prior can be used to set a lower and upper bound for the mean proportion of
#'   partitions lost. Note that for proportions close to 0 or to
#'   `partition_loss_max` (see below), the resulting mean partition loss can
#'   slightly differ from what is specified here, because we internally
#'   translate this prior to the logit scale.
#' @param partition_loss_mean_prior_upper Prior (95% quantile) for the mean
#'   relative partition loss (see `partition_loss_mean_prior_lower` for
#'   details). In well-functioning dPCR assays, the mean proportion of lost
#'   partitions should not be very high (definitely below 50%).
#' @param partition_loss_variation_prior_lower Prior (5% quantile) for the
#'   variation in the number of invalid partitions across dPCR runs. The
#'   proportion of partitions lost typically varies between dPCR runs. We thus
#'   model this proportion as logit-normal distributed, with mean
#'   (approximately) defined by `partition_loss_mean_prior` and logit-level
#'   standard deviation `sigma`. You can use
#'   `partition_loss_variation_prior_lower` and
#'   `partition_loss_variation_prior_upper` to set a lower and upper bound for
#'   `sigma`.
#' @param partition_loss_variation_prior_upper Prior (95% quantile) for the
#'   variation in the number of invalid partitions across dPCR runs (see
#'   `partition_loss_variation_prior_lower` for details).
#' @param partition_loss_max The maximum proportion of partitions that can be
#'   lost in a valid dPCR run. During quality control, runs where the proportion
#'   of invalid partitions is above some threshold (e.g. 50%) are often
#'   discarded. This parameter can be used to represent such a QC threshold.
#' @param volume_scaled_prior_mu Prior (mean) on the conversion factor
#'   (partition volume scaled by the dilution of wastewater in the assay) for
#'   the dPCR reaction. See details for further explanation.
#' @param volume_scaled_prior_sigma Prior (standard deviation) on the conversion
#'   factor (partition volume scaled by the dilution of wastewater in the assay)
#'   for the dPCR reaction. If this is set to zero, the conversion factor will
#'   be fixed to the prior mean and not estimated.
#' @param pre_replicate_cv_prior_mu Prior (mean) on the coefficient of variation
#'   of concentrations *before* the replication stage.
#' @param pre_replicate_cv_prior_sigma Prior (standard deviation) on the
#'   coefficient of variation of concentrations *before* the replication stage.
#' @param prePCR_noise_type The parametric distribution to assume for noise
#'   before the PCR assay. Currently supported are "log-normal" and "gamma". The
#'   choice of the parametric distribution typically makes no relevant
#'   difference for the noise model, but can make a relevant difference for the
#'   LOD model if [LOD_estimate_dPCR()] is used.
#'
#' @param use_taylor_approx If TRUE (default), a Taylor expansion approximation
#'   is used to estimate the CV of measurements under pre-PCR noise. The
#'   approximation is very accurate, unless concentrations are extremely high
#'   (so high that the quality of the measurements from dPCR would anyway be
#'   questionable).
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @keywords internal
noise_estimate_ <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           cv_type = "constant",
           total_partitions_observe = NULL,
           max_partitions_prior_lower= NULL,
           max_partitions_prior_upper = NULL,
           partition_loss_mean_prior_lower = NULL,
           partition_loss_mean_prior_upper = NULL,
           partition_loss_variation_prior_lower = NULL,
           partition_loss_variation_prior_upper = NULL,
           partition_loss_max = NULL,
           volume_scaled_prior_mu = NULL,
           volume_scaled_prior_sigma = NULL,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           prePCR_noise_type = "log-normal",
           use_taylor_approx = TRUE,
           modeldata = modeldata_init()) {

    if (!is.null(modeldata$obs_dist) && (modeldata$obs_dist == 4 && cv_type != "dPCR")) {
      cli::cli_abort(paste0(
        "You specified positive partitions from a dPCR assay as measurements, ",
        "but the noise model is not compatible with this observation type. ",
        "Please use `noise_estimate_dPCR()` instead."
      ))
    }

    modeldata$pr_noise <- replicates

    modeldata$nu_upsilon_a_prior <- set_prior(
      "nu_upsilon_a", "truncated normal",
      mu = cv_prior_mu, sigma = cv_prior_sigma
    )
    modeldata$.init$nu_upsilon_a <- 0.1 # 10% coefficient of variation

    if (cv_type == "constant") {
      modeldata$total_partitions_observe <- FALSE
      modeldata$cv_type <- 0
      modeldata$max_partitions_prior <- numeric(0)
      modeldata$partition_loss_mu_prior <- numeric(0)
      modeldata$partition_loss_sigma_prior <- numeric(0)
      modeldata$partition_loss_max <- numeric(0)
      modeldata$.init$max_partitions <- numeric(0)
      modeldata$.init$partition_loss_mu <- numeric(0)
      modeldata$.init$partition_loss_sigma <- numeric(0)
      modeldata$.init$partition_loss_raw <- numeric(0)
      modeldata$nu_upsilon_c_prior <- numeric(0)
      modeldata$.init$nu_upsilon_c <- numeric(0)
      modeldata$cv_pre_type <- numeric(0)
      modeldata$cv_pre_approx_taylor <- numeric(0)
    } else if (cv_type == "dPCR") {
      # 1: continuous model of measured concentrations
      # 3: binomial model of positive partitions
      modeldata$cv_type <- tbe(
        ifelse(modeldata$obs_dist == 4, 3, 1), "obs_dist"
      )

      if (total_partitions_observe || (!is.null(modeldata$obs_dist) && modeldata$obs_dist == 4)) {
        modeldata$total_partitions_observe <- TRUE
        modeldata$max_partitions_prior <- numeric(0)
        modeldata$partition_loss_mu_prior <- numeric(0)
        modeldata$partition_loss_max <- numeric(0)
        modeldata$partition_loss_sigma_prior <- numeric(0)
        modeldata$.init$max_partitions <- numeric(0)
        modeldata$.init$partition_loss_mu <- numeric(0)
        modeldata$.init$partition_loss_sigma <- numeric(0)
        modeldata$.init$partition_loss_raw <- numeric(0)
        modeldata$.checks$check_total_partitions_col <- function(md, ...) {
          if (length(md$dPCR_total_partitions)==0) {
            cli::cli_abort(paste0(
              "You specified `total_partitions_observe = TRUE`, which requires ",
              "a column with the number of total partitions in the PCR for ",
              "each sample in your data. Please specify such a column via the ",
              "`total_partitions_col` argument in ",
              cli_help("concentrations_observe"), "."
            ))
          }
        }
      } else {
        modeldata$total_partitions_observe <- FALSE

        # maximum number of partitions
        modeldata$max_partitions_prior <- set_prior_trunc_normal(
          "max_partitions", "truncated normal",
          q5 = max_partitions_prior_lower * 1e-4, # scale by 1e-4 for numerical efficiency
          q95 = max_partitions_prior_upper * 1e-4
        )
        modeldata$.init$max_partitions <- init_from_location_scale_prior(
          modeldata$max_partitions_prior, enforce_positive = TRUE
        )

        # mean partition loss
        modeldata$partition_loss_mu_prior <- set_prior_normal(
          "partition_loss_mu", "truncated normal",
          q5 = qlogis(partition_loss_mean_prior_lower/partition_loss_max),
          q95 = qlogis(partition_loss_mean_prior_upper/partition_loss_max)
        )
        modeldata$.init$partition_loss_mu <- init_from_location_scale_prior(
          modeldata$partition_loss_mu_prior
        )

        # variation in partition loss
        modeldata$partition_loss_sigma_prior <- set_prior_trunc_normal(
          "partition_loss_sigma", "truncated normal",
          q5 = partition_loss_variation_prior_lower,
          q95 = partition_loss_variation_prior_upper
        )
        modeldata$.init$partition_loss_sigma <- init_from_location_scale_prior(
          modeldata$partition_loss_sigma_prior, enforce_positive = TRUE
        )

        # maximum partition loss (threshold)
        modeldata$partition_loss_max <- partition_loss_max

        # non-centered noise for partition loss
        modeldata$.init$partition_loss_raw <- tbe(
          rep(-1e-4, sum(modeldata$n_averaged)),
          "n_averaged"
        )

      }

      # conversion factor for dPCR
      modeldata$nu_upsilon_c_prior <- set_prior(
        "nu_upsilon_c", "truncated normal",
        mu = volume_scaled_prior_mu * 1e+5, # scaling by 1e+5 for numerical reasons
        sigma = volume_scaled_prior_sigma * 1e+5
      )
      modeldata$.init$nu_upsilon_c <- init_from_location_scale_prior(
        modeldata$nu_upsilon_c_prior, enforce_positive = TRUE
        )

      if (prePCR_noise_type == "gamma") {
        modeldata$cv_pre_type <- 0
      } else if (prePCR_noise_type %in% c("log-normal", "lognormal")) {
        modeldata$cv_pre_type <- 1
      } else {
        cli::cli_abort(paste0(
            "`prePCR_noise_type = ", prePCR_noise_type, "` not supported.",
            "Available options: 'gamma', `log-normal`."
          ))
      }
      modeldata$cv_pre_approx_taylor <- use_taylor_approx

    } else if (cv_type == "constant_var") {
      modeldata$total_partitions_observe <- FALSE
      modeldata$cv_type <- 2
      modeldata$max_partitions_prior <- numeric(0)
      modeldata$partition_loss_mu_prior <- numeric(0)
      modeldata$partition_loss_sigma_prior <- numeric(0)
      modeldata$partition_loss_max <- numeric(0)
      modeldata$.init$max_partitions <- numeric(0)
      modeldata$.init$partition_loss_mu <- numeric(0)
      modeldata$.init$partition_loss_sigma <- numeric(0)
      modeldata$.init$partition_loss_raw <- numeric(0)
      modeldata$nu_upsilon_c_prior <- numeric(0)
      modeldata$.init$nu_upsilon_c <- numeric(0)
      modeldata$cv_pre_type <- numeric(0)
      modeldata$cv_pre_approx_taylor <- numeric(0)
    } else {
      cli::cli_abort(
        paste0(
          "Noise type `", cv_type, "` not supported. Available options: ",
          "'constant', `dPCR`, `constant_var`."
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

      modeldata$.checks$check_replicate_ids <- function(md, ...) {
        if (!"replicate_ids" %in% names(md)) {
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

    if (cv_type == "constant") {
      modeldata$.str$measurements[["noise"]] <- list(
        noise_estimate = .str_details
      )
    } else if (cv_type == "dPCR") {
      modeldata$.str$measurements[["noise"]] <- list(
        noise_estimate_dPCR = .str_details
      )
    } else if (cv_type == "constant_var") {
      modeldata$.str$measurements[["noise"]] <- list(
        noise_estimate_constant_var = .str_details
      )
    }

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
#'   [noise_estimate_dPCR()].
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
#' @family {noise models}
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

#' Estimate measurement noise for digital PCR data
#'
#' @description This option estimates the unexplained variation in wastewater
#'   measurements using a coefficient of variation model specialized for digital
#'   PCR (e.g. digital droplet PCR). Specifically, the coefficient of variation
#'   is modeled as a function of the expected concentration according to the
#'   statistical properties of dPCR.
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
#'   the total number of partitions and conversion factor). Moreover, when
#'   `replicates=TRUE`, this is only the CV after the replication stage (see
#'   details for more explanation).
#'
#' @details The conversion factor (see `volume_scaled_prior_mu`,
#'   `volume_scaled_prior_sigma`) is the partition volume v multiplied with a
#'   scaling factor s. The scaling factor accounts for concentration differences
#'   between the sample and the reaction mix, for example due to extraction or
#'   adding of reagents. For example, if the partition volume is 4.5e-7 mL and
#'   the scaling factor is 100:3 (i.e. 100 gc/mL in the original sample
#'   correspond to 3 gc/mL in the PCR reaction), then the overall conversion
#'   factor is 4.5e-7 * 100 / 3 = 1.5e-5.
#'
#' @details When `replicates=TRUE`, two coefficients of variation are estimated:
#' - the CV before the replication stage (see `pre_replicate_cv_prior_mu`)
#' - the CV after the replication stage (see `cv_prior_mu`)
#'
#'   The meaning of these CV estimates depends on the type of replicates. If the
#'   replicates are biological replicates (i.e. independently processed), then
#'   `cv` estimates the noise in the preprocessing before the PCR, and
#'   `pre_replicate_cv` estimates the noise from anything before preprocessing
#'   (e.g. sampling noise and all other unexplained variation). In contrast, if
#'   the replicates are technical replicates (i.e. several PCR runs of the same
#'   preprocessed sample), then `cv` estimates only unexplained PCR noise
#'   (should be close to zero), and `pre_replicate_cv` estimates all other noise
#'   (including preprocessing noise.)
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements (`cv`): `Truncated normal`
#' - maximum number of total partitions in dPCR: `Truncated normal`
#' - mean proportion of lost partitions in dPCR: `Normal (logit-level)`
#' - variation of proportion of lost partitions: `Truncated normal (logit-level)`
#' - conversion factor for dPCR: `Truncated normal`
#' - coefficient of variation of concentration before the replication stage (`pre_replicate_cv`):
#'   `Truncated normal`
#'
#' @inheritParams noise_estimate_
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {noise models}
#' @seealso [LOD_estimate_dPCR] for a limit of detection model specialised for
#'   dPCR.
noise_estimate_dPCR <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           total_partitions_observe = FALSE,
           max_partitions_prior_lower = 5000,
           max_partitions_prior_upper = 30000,
           partition_loss_mean_prior_lower = 0.01,
           partition_loss_mean_prior_upper = 0.3,
           partition_loss_variation_prior_lower = 0.5,
           partition_loss_variation_prior_upper = 2,
           partition_loss_max = 0.5,
           volume_scaled_prior_mu = 1e-5,
           volume_scaled_prior_sigma = 4e-5,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           prePCR_noise_type = "log-normal",
           use_taylor_approx = TRUE,
           modeldata = modeldata_init()) {
    return(noise_estimate_(
      replicates = replicates,
      cv_prior_mu = cv_prior_mu,
      cv_prior_sigma = cv_prior_sigma,
      cv_type = "dPCR",
      max_partitions_prior_lower = max_partitions_prior_lower,
      max_partitions_prior_upper = max_partitions_prior_upper,
      partition_loss_mean_prior_lower = partition_loss_mean_prior_lower,
      partition_loss_mean_prior_upper = partition_loss_mean_prior_upper,
      partition_loss_variation_prior_lower = partition_loss_variation_prior_lower,
      partition_loss_variation_prior_upper = partition_loss_variation_prior_upper,
      partition_loss_max = partition_loss_max,
      total_partitions_observe = total_partitions_observe,
      volume_scaled_prior_mu = volume_scaled_prior_mu,
      volume_scaled_prior_sigma = volume_scaled_prior_sigma,
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
#'   see [noise_estimate_dPCR()].
#'
#' @description If multiple measurements (replicates) per sample are provided,
#'   `EpiSewer` can also explicitly model variation before the replication
#'   stage.
#'
#' @details Note that although this model keeps the variance constant, the prior
#'   for the measurement noise is still in terms of the (average) coefficient of
#'   variation (CV). This makes prior specification easier since the CV is
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
#' @family {noise models}
noise_estimate_constant_var <-
  function(replicates = FALSE,
           cv_prior_mu = 0,
           cv_prior_sigma = 1,
           pre_replicate_cv_prior_mu = 0,
           pre_replicate_cv_prior_sigma = 1,
           warn = TRUE,
           modeldata = modeldata_init()) {
    if (warn) {
      cli::cli_inform(c("!" = paste0(
      "You have specified ",
      cli_help("noise_estimate_constant_var"),
      " as the model component for measurement noise.",
      " Note that modeling a constant variance",
      " is likely a model misspecification and should only be used for",
      " comparison purposes with better models like ",
      cli_help("noise_estimate"), " or ", cli_help("noise_estimate_dPCR"), ".",
      " You can specify ",
      "{.code noise_estimate_constant_var(warn=FALSE)} to disable this warning."
      )))
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
  if (!is.null(modeldata$obs_dist) && modeldata$obs_dist == 4) {
    cli::cli_abort(paste0(
      "You specified positive partitions from a dPCR assay as measurements, ",
      "which means that non-detects (zero positive partitions) must be ",
      "modeled. Please use `LOD = LOD_dPCR()`."
    ))
  }

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
#'   exponential model can be derived from the statistical properties of dPCR,
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
#' @details If used together with [noise_estimate_dPCR()], EpiSewer will also
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

  if (!LOD_type %in% c("exponential", "dPCR")) { # "dPCR" is synonym
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

#' Estimate a limit of detection model for digital PCR data
#'
#' @description Pathogen concentrations below a certain threshold may not be
#'   detectable and thus erroneously measured as 0. This option adjusts for a
#'   limit of detection based on the statistical properties of digital PCR
#'   (dPCR) and includes zero measurements in the likelihood.
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
#'   uses the number of partitions in the dPCR reaction and the conversion
#'   factor as defined and estimated by [noise_estimate_dPCR()]. It can
#'   therefore only be used together with `noise = noise_estimate_dPCR()` in
#'   [model_measurements()].
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @family {LOD models}
LOD_estimate_dPCR <- function(drop_prob = 1e-10, modeldata = modeldata_init()) {

  modeldata <- tbc("LOD_estimate_dPCR",
    {
      if (is.null(modeldata$cv_type) || !(modeldata$cv_type %in% c(1,3))) {
        cli::cli_abort(paste0(
          "To use LOD = ",
          cli_help("LOD_estimate_dPCR"), ", you must specify noise = ",
          cli_help("noise_estimate_dPCR"), ' in ',
          cli_help("model_measurements"), "."
        ))
      }

      if (!is.null(modeldata$obs_dist) && modeldata$obs_dist == 4) {
        # when using a binomial model of positive partitions,
        # non-detects are automatically accounted for and don't have to be modeled
        modeldata$LOD_model <- 0
        modeldata$LOD_scale <- numeric(0)
        modeldata$LOD_drop_prob <- 0
        return(modeldata)
      } else {
        modeldata$LOD_model <- 2
        modeldata$LOD_scale <- numeric(0)
        modeldata$LOD_drop_prob <- drop_prob
        return(modeldata)
      }
    },
    required = c(
      "cv_type"
      ),
    modeldata = modeldata
    )

  modeldata$LOD_drop_prob <- drop_prob

  modeldata$.str$measurements[["LOD"]] <- list(
    LOD_estimate_dPCR = c()
  )

  return(modeldata)
}
