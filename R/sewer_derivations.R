# Registry of derived modeldata variables.
#
# Cross-cutting variables that depend on the outputs of several model
# components are not computed inside component bodies but registered here as
# small pure derivation units. They are resolved on demand via md_need()
# (with memoization into the modeldata and cycle detection) and by the settle
# pass at the end of modeldata_compile(). This replaces the v1
# modeldata_update_metainfo() fixpoint machinery.

#' Create a derivation unit
#'
#' @param provides Character vector of `$`-paths this unit provides.
#' @param requires Character vector of `$`-paths this unit needs (resolved
#'   recursively via [md_need()] before the unit runs). Arguments of model
#'   components (available via [md_spec()]) are always available and need not
#'   be declared.
#' @param fn `function(md)` returning a named list keyed by the `provides`
#'   paths. Must be pure apart from cli messages (which fire exactly once,
#'   thanks to memoization).
#' @keywords internal
md_derive <- function(provides, requires, fn) {
  list(provides = provides, requires = requires, fn = fn)
}

#' The registry of all derivation units
#'
#' @description Returns the list of derivation units. Each derived variable
#'   is provided by exactly one unit; [md_need()] resolves them on demand.
#' @keywords internal
modeldata_derivations <- function() {
  list(
    length_shedding = md_derive(
      provides = ".metainfo$length_shedding",
      requires = c("S", "D", "T"),
      fn = function(md) {
        list(".metainfo$length_shedding" = md$S + md$D + md$T)
      }
    ),

    length_I = md_derive(
      provides = ".metainfo$length_I",
      requires = c("L", "S", "D", "T"),
      fn = function(md) {
        list(".metainfo$length_I" = md$L + md$S + md$D + md$T)
      }
    ),

    total_delay_dist = md_derive(
      provides = ".metainfo$total_delay_dist",
      requires = c("residence_dist", "shedding_dist", "incubation_dist"),
      fn = function(md) {
        inc_shed_dist <- convolve(
          md$incubation_dist, rev(md$shedding_dist), type = "o"
        )
        total_delay_dist <- convolve(
          inc_shed_dist, rev(md$residence_dist), type = "o"
        )
        list(".metainfo$total_delay_dist" = total_delay_dist)
      }
    ),

    partial_window = md_derive(
      provides = ".metainfo$partial_window",
      requires = c(".metainfo$total_delay_dist"),
      fn = function(md) {
        list(".metainfo$partial_window" = which(
          cumsum(md$.metainfo$total_delay_dist) > 0.9
        )[1] - 1) # delay dist is 0-indexed
      }
    ),

    partial_generation = md_derive(
      provides = ".metainfo$partial_generation",
      requires = c("generation_dist"),
      fn = function(md) {
        list(".metainfo$partial_generation" = which(
          cumsum(md$generation_dist) > 0.9
        )[1]) # generation dist is 1-indexed
      }
    ),

    LOD_expected_scale = md_derive(
      provides = ".metainfo$LOD_expected_scale",
      requires = c(
        "LOD_model", "LOD_scale", "n_averaged",
        "dPCR_total_partitions", "total_partitions_observe",
        "max_partitions_prior", "partition_loss_mu_prior",
        "nu_upsilon_c_prior"
      ),
      fn = function(md) {
        if (md$LOD_model == 0) {
          scale <- NA
        } else if (md$LOD_model == 1) {
          scale <- md$LOD_scale
        } else if (md$LOD_model == 2) {
          total_partitions_median <- median(md$dPCR_total_partitions)
          n_averaged_median <- median(md$n_averaged)

          if (md$total_partitions_observe) {
            total_partitions_expected <- total_partitions_median
          } else {
            max_partitions_expected <- 1 / 2 * (
              md$max_partitions_prior$max_partitions_prior[1] +
                md$max_partitions_prior$max_partitions_prior[2]
            )
            total_partitions_expected <- max_partitions_expected *
              (1 - plogis(
                md$partition_loss_mu_prior$partition_loss_mu_prior[1]
              ))
          }

          conversion_expected <- 1 / 2 * (
            md$nu_upsilon_c_prior$nu_upsilon_c_prior[1] +
              md$nu_upsilon_c_prior$nu_upsilon_c_prior[2]
          )

          scale <- (
            total_partitions_expected * conversion_expected * n_averaged_median
          )
        } else {
          scale <- NA
        }
        list(".metainfo$LOD_expected_scale" = scale)
      }
    ),

    initial_cases_crude = md_derive(
      provides = ".metainfo$initial_cases_crude",
      requires = c(
        "measured_concentrations", "measure_to_sample", "sample_to_date",
        "flow", ".metainfo$composite_window", ".metainfo$load_per_case"
      ),
      fn = function(md) {
        # crude descriptive estimate of cases at start of time series
        # here we take the mean of the first week of samples
        initial_cases_crude <-
          0.1 + # small offset to avoid zero cases
          mean(
            md$measured_concentrations[
              md$measure_to_sample[which(md$sample_to_date <= 7)]
            ],
            na.rm = TRUE
          ) *
            mean(md$flow[1:7], na.rm = TRUE) /
            md$.metainfo$load_per_case
        list(".metainfo$initial_cases_crude" = initial_cases_crude)
      }
    ),

    load_curve_crude = md_derive(
      provides = ".metainfo$load_curve_crude",
      requires = c(
        "measured_concentrations", "measure_to_sample", "sample_to_date",
        "flow", ".metainfo$total_delay_dist", ".metainfo$length_I", "T",
        ".metainfo$T_start_date", ".metainfo$LOD_expected_scale"
      ),
      fn = function(md) {
        list(".metainfo$load_curve_crude" = get_load_curve_crude(
          md$measured_concentrations, md$measure_to_sample,
          md$sample_to_date, md$flow,
          md$.metainfo$total_delay_dist,
          max_shift = md$.metainfo$length_I - md$T,
          md$.metainfo$T_start_date,
          # asymptotic posterior expectation for non-detects
          impute_zero = 1 / md$.metainfo$LOD_expected_scale,
          impute_zero_runs = TRUE,
          interpolate = TRUE, loess_window = 56,
          plot_smoothed_curve = FALSE
        ))
      }
    ),

    infection_curve_crude = md_derive(
      provides = ".metainfo$infection_curve_crude",
      requires = c(".metainfo$load_curve_crude", ".metainfo$load_per_case"),
      fn = function(md) {
        icc <- get_infection_curve_crude(
          md$.metainfo$load_curve_crude, md$.metainfo$load_per_case
        )
        if (max(icc$infections, na.rm = TRUE) < 1) {
          cli::cli_abort(c("!" = paste(
            "With the current data and model specification, the estimated number",
            "of infections will be too low",
            paste0(
              "(maximum daily incidence: ",
              round(max(icc$infections, na.rm = TRUE), 4),
              ")."
            ),
            "Please ensure that the volume units of the concentrations and flows",
            "are identical (ideally in mL) and check the assumed `load_per_case`."
          )))
        } else if (max(icc$infections, na.rm = TRUE) < 5) {
          cli::cli_inform(c("!" = paste(
            "With the current data and model specification, the estimated number",
            "of infections will be quite low",
            paste0(
              "(maximum daily incidence: ",
              round(max(icc$infections, na.rm = TRUE), 4),
              ")."
            ),
            "Please ensure that the volume units of the concentrations and flows",
            "are identical (ideally in mL) and check the assumed `load_per_case`."
          )))
        } else if (median(icc$infections, na.rm = TRUE) > 500000) {
          cli::cli_inform(c("!" = paste(
            "With the current data and model specification, the estimated number",
            "of infections will be extremely high",
            paste0(
              "(median daily incidence: ",
              round(median(icc$infections, na.rm = TRUE), 0),
              ")."
            ),
            "Please ensure that the volume units of the concentrations and flows",
            "are identical (ideally in mL) and check the assumed `load_per_case`."
          )))
        }
        list(".metainfo$infection_curve_crude" = icc)
      }
    ),

    # Seeding extension: determines the length of the seeding phase from
    # non-detects at the start of the measurement time series, and with it
    # the lengths of the modeled Rt time series. Provides the stan variable
    # `se` (seeding extension in days).
    seeding_extension = md_derive(
      provides = c(
        "se", ".metainfo$date_triple_detect", ".metainfo$length_seeding",
        ".metainfo$length_R", ".metainfo$length_R_modeled"
      ),
      requires = c(
        ".metainfo$load_curve_crude", "L", "S", "D", "T", "G",
        ".metainfo$T_start_date"
      ),
      fn = function(md) {
        extend_seeding <- md_spec(md, "seeding", "extend", default = TRUE)
        seeding_helper <- md_spec(
          md, "seeding", ".helper", default = "seeding_estimate_rw"
        )

        triplets <- md$.metainfo$load_curve_crude[detect_next_n >= 3]
        if (nrow(triplets) == 0) {
          date_triple_detect <- NA
          seed_extension <- 0
          cli::cli_inform(c("!" = paste(
            "The measurement time series contains a large percentage of",
            "non-detects (zero concentration measurements).",
            "EpiSewer will attempt to model the transmission dynamics,",
            "but there could be sampling problems.",
            "Please make sure to check the estimated infection time series",
            "after model fitting, as infection incidence could be very low."
          )))
        } else {
          date_triple_detect <- triplets[, min(date)]
          seed_length <- as.numeric(
            date_triple_detect - md$.metainfo$T_start_date
          )
          if (seed_length > md$G && extend_seeding) {
            seed_extension <- seed_length - md$G
          } else {
            seed_extension <- 0
          }
        }

        if (seed_extension > 0) {
          cli::cli_inform(c("i" = paste0(
            "Due to non-detects at the start of the measurement time series, ",
            "Rt will only be estimated from ",
            date_triple_detect, " onwards. ",
            "Measurements and Rt before that date are still modeled, but using ",
            "an extended seeding phase. Set 'seeding = ",
            seeding_helper, "(extend = FALSE)' ",
            "to deactivate this behaviour."
          )))
        }

        length_R <- md$L + md$S + md$D + md$T - md$G

        list(
          "se" = seed_extension,
          ".metainfo$date_triple_detect" = date_triple_detect,
          ".metainfo$length_seeding" = md$G + seed_extension,
          ".metainfo$length_R" = length_R,
          ".metainfo$length_R_modeled" = length_R - seed_extension
        )
      }
    ),

    # Reconciled decision whether the total number of dPCR partitions is
    # observed from the data. Rooted purely in component arguments: positive
    # partition observations force it to TRUE, otherwise the noise component's
    # argument decides (FALSE for non-dPCR noise models, which lack the
    # argument).
    total_partitions_observe = md_derive(
      provides = "total_partitions_observe",
      requires = character(0),
      fn = function(md) {
        observation_type <- md_spec(
          md, "concentrations", "observation_type", default = "concentrations"
        )
        if (identical(observation_type, "partitions")) {
          value <- TRUE
        } else {
          value <- isTRUE(md_spec(
            md, "noise", "total_partitions_observe", default = FALSE
          ))
        }
        list("total_partitions_observe" = value)
      }
    )
  )
}
