#' Model the shedding process
#'
#' @description This module function is used to specify the components of the
#'  `shedding` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param incubation_dist Incubation period distribution. `EpiSewer` uses this
#'   as a proxy for the time between infection and the start of shedding, as
#'   shedding load distributions in the literature are often given from symptom
#'   onset onwards. If the assumed shedding load distribution instead starts
#'   from the time of infection, the incubation period should be fixed to 0
#'   days. Modeling options:
#' `r component_functions_("incubation_dist")`
#' @param shedding_dist Shedding load distribution. Describes how the total load
#'   shed by an individual is distributed over time (and therefore sums to 1).
#'   Modeling options:
#' `r component_functions_("shedding_dist")`
#' @param load_per_case Average total load per case. This is a scaling factor
#'   that describes how many pathogen particles are shed by the average infected
#'   individual overall and how much of this is detectable at the sampling site.
#'   It depends both on biological factors as well as on the specific
#'   sewage system. Modeling options:
#' `r component_functions_("load_per_case")`
#' @param load_variation Individual-level shedding load variation. The strength
#' of shedding may vary between individuals. Modeling this variation can better
#' account for uncertainty especially at low incidence. Modeling options:
#' `r component_functions_("load_variation")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `shedding` module.
#' @export
#' @family {module functions}
model_shedding <- function(
    incubation_dist = incubation_dist_assume(),
    shedding_dist = shedding_dist_assume(),
    load_per_case = load_per_case_calibrate(),
    load_variation = load_variation_estimate()) {
  verify_is_modeldata(incubation_dist, "incubation_dist")
  verify_is_modeldata(shedding_dist, "shedding_dist")
  verify_is_modeldata(load_per_case, "load_per_case")
  verify_is_modeldata(load_variation, "load_variation")
  modeldata <- modeldata_combine(
    incubation_dist, shedding_dist, load_per_case, load_variation
  )
  return(modeldata)
}

#' Assume an incubation period distribution
#'
#' @description This option assumes a fixed incubation period distribution for
#'   the shedding model in `EpiSewer`.
#'
#' @details The incubation period is the time between infection and symptom
#'   onset. This assumption is used when `shedding_reference="symptom_onset"`,
#'   i.e. to support shedding load distributions referenced by days since
#'   symptom onset.
#'
#' @param incubation_dist A numeric vector representing a discrete incubation
#'   period distribution, starting with the probability for an incubation period
#'   of 0 days, 1 day, 2 days, and so on.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso Helpers to discretize continuous probability distributions:
#'   [get_discrete_gamma()], [get_discrete_lognormal()]
incubation_dist_assume <-
  function(incubation_dist = NULL, modeldata = modeldata_init()) {

    modeldata <- tbc("incubation_dist_compute", {
      if (modeldata$.metainfo$shedding_reference == "symptom_onset") {

        modeldata <- tbp("incubation_dist_assume",
         {
           modeldata$L <- length(incubation_dist) - 1
           incubation_dist <- check_dist(
             incubation_dist, "incubation period distribution"
           )
           modeldata$incubation_dist <- incubation_dist
           return(modeldata)
         },
         required_assumptions = "incubation_dist",
         modeldata = modeldata
        )

      } else if (modeldata$.metainfo$shedding_reference == "infection") {
        # this is a workaround because the incubation period is currently
        # not needed for anything else than for modeling the shedding
        # profile.
        modeldata$L <- 1
        modeldata$incubation_dist <- c(1, 0)
        return(modeldata)
      }

    },
    required = c(
      ".metainfo$shedding_reference"
    ),
    modeldata = modeldata
    )

    modeldata$.str$shedding[["incubation_dist"]] <- list(
      incubation_dist_assume = c()
    )

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
#' @param shedding_reference Is the shedding load distribution relative to the
#'   day of `"infection"` or the day of `"symptom_onset"`? This is important because
#'   shedding load distributions provided in the literature are sometimes by
#'   days since infection and sometimes by days since symptom onset. If
#'   `shedding_reference="symptom_onset"`, EpiSewer also needs information about
#'   the incubation period distribution (see [incubation_dist_assume()]).
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso Helpers to discretize continuous probability distributions:
#'   [get_discrete_gamma()], [get_discrete_lognormal()]
shedding_dist_assume <-
  function(shedding_dist = NULL, shedding_reference = NULL, modeldata = modeldata_init()) {
    modeldata <- tbp("shedding_dist_assume",
      {
        modeldata$S <- length(shedding_dist) - 1
        shedding_dist <- check_dist(shedding_dist, "shedding load distribution")
        modeldata$shedding_dist <- shedding_dist
        if (!shedding_reference %in% c("infection", "symptom_onset")) {
          cli::cli_abort(paste(
            "The provided `shedding_reference` argument is invalid.",
            'Must be either "infection" or "symptom_onset".'
          ))
        }
        modeldata$.metainfo$shedding_reference <- shedding_reference

        modeldata$.str$shedding[["shedding_dist"]] <- list(
          shedding_dist_assume = c(shedding_reference = shedding_reference)
        )

        return(modeldata)
      },
      required_assumptions = c("shedding_dist", "shedding_reference"),
      modeldata = modeldata
    )

    modeldata$.str$shedding[["shedding_dist"]] <- list(
      shedding_dist_assume = c()
    )

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
#' @seealso {Helper for finding a suitable load per case assumption:}
#'   [suggest_load_per_case()]
#' @family {load per case functions}
load_per_case_assume <-
  function(load_per_case = NULL, modeldata = modeldata_init()) {
    modeldata <- tbp("load_per_case_assume",
      {
        modeldata$load_mean <- load_per_case
        modeldata$.metainfo$load_per_case <- load_per_case
        return(modeldata)
      },
      required_assumptions = "load_per_case",
      modeldata = modeldata
    )

    modeldata$.str$shedding[["load_per_case"]] <- list(
      load_per_case_assume = c()
    )

    return(modeldata)
  }

#' Calibrate the average load per case using case count data
#'
#' @description This option calibrates the average total shedding load per case
#'   based on the relationship between measured concentrations and case counts.
#'   The goal is to obtain a `load_per_case` that is on the right order of
#'   magnitude - this will not be sufficient for accurate prevalence estimation
#'   from wastewater, but fully sufficient for monitoring trends and estimating
#'   Rt.
#'
#' @param cases A `data.frame` of case numbers with each row representing one
#'   day. Must have at least a column with dates and a column with case numbers.
#' @param min_cases This is an alternative to supplying a `data.frame` of cases.
#'   If `min_cases` is specified, then `load_per_case` is calibrated based on
#'   the smallest observed load in the wastewater. EpiSewer uses `min_cases =
#'   10` as a default if no case data is supplied (see [sewer_assumptions()]).
#' @param ascertainment_prop Proportion of all cases that get detected /
#'   reported. Can be used to account for underreporting of infections. Default
#'   is `ascertainment_prop=1`, meaning that 100% of infections become confirmed
#'   cases.
#' @param measurement_shift The specific timing between wastewater
#'   concentrations and case numbers depends on reporting delays and shedding
#'   profiles and is typically uncertain. This argument allows to shift the
#'   concentration and case number time series relative to each other and to
#'   average over several potential lags/leads, as specified by an integer
#'   vector. The default is `measurement_shift = seq(-7,7)`, i.e. a shift of
#'   concentrations between up to one week before and after case numbers.
#' @param shift_weights Weights for the shifted comparisons. Must be an numeric
#'   vector of the same length as `measurement_shift`. If `NULL` (default), the
#'   weights are chosen to be approximately inversely proportional to the shift
#'   distance.
#' @param date_col Name of the date column in the provided cases `data.frame`.
#' @param case_col Name of the column containing the case numbers.
#' @param signif_fig Significant figures to round to. Since this heuristic only
#'   provides crude estimates which should not be over-interpreted, the result
#'   gets rounded (this also increases stability when fitting the model at
#'   different points in time). Default is rounding to the 2 most significant
#'   figures.
#'
#' @details In `EpiSewer`, the `load_per_case` serves as a scaling factor
#'   describing how many pathogen particles are shed by the average infected
#'   individual overall and how much of this is detectable at the sampling site.
#'   This depends both on biological factors as well as on the specific sewage
#'   system and laboratory quantification. It is therefore almost always
#'   necessary to assume the load per case based on a comparison of measured
#'   concentrations/loads and case numbers.
#'
#' @details If a `data.frame` of cases is supplied via the `cases` argument, the
#'   average load per case is determined by fitting a linear regression model
#'   with loads (computed using concentrations and flows) as dependent variable
#'   and case numbers as independent variable. This does not explicitly account
#'   for shedding profiles or reporting delays, but the `measurement_shift`
#'   argument allows to average over a set of relative shifts between the two
#'   time series.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {load per case functions}
load_per_case_calibrate <- function(cases = NULL, min_cases = NULL,
                                    ascertainment_prop = 1,
                                    measurement_shift = seq(-7,7),
                                    shift_weights = 1/(abs(measurement_shift)+1),
                                    date_col = "date",
                                    case_col = "cases",
                                    signif_fig = 2,
                                    modeldata = modeldata_init()) {
  modeldata <- tbp("load_per_case_calibrate",
   {
     if (!is.null(cases)) {
       modeldata <- tbc("load_per_case_cases_suggest", {
         measurements <- setnames(
           copy(modeldata$.sewer_data$concentrations_observe$measurements),
           old = with(
             modeldata$.metainfo$measurements_cols,
             c(date_col, concentration_col)
           ),
           new = c("date", "concentration")
         )

         flows <- setnames(
           copy(modeldata$.sewer_data$flows_observe$flows),
           old = with(
             modeldata$.metainfo$flows_cols,
             c(date_col, flow_col)
           ),
           new = c("date", "flow")
         )

         required_data_cols <- c(date_col, case_col)
         if (!all(required_data_cols %in% names(cases))) {
           cli::cli_abort(
             c(paste(
               "The following columns must be present",
               "in the provided cases `data.frame`:",
               paste(required_data_cols, collapse = ", ")
             ),
             paste("Please adjust the `data.frame` or specify the right column",
                   "names via the `_col` arguments of this function."))
           )
         }
         cases = as.data.table(cases)[, .SD, .SDcols = required_data_cols]
         setnames(cases, old = required_data_cols, new = c("date", "cases"))

         suggested_load <- suggest_load_per_case(
           measurements = measurements,
           cases = cases,
           flows = flows,
           ascertainment_prop = ascertainment_prop,
           measurement_shift = measurement_shift,
           shift_weights = shift_weights,
           date_col = "date",
           concentration_col = "concentration",
           flow_col = "flow",
           case_col = "cases",
           signif_fig = signif_fig
         )
         modeldata$load_mean <- suggested_load
         modeldata$.metainfo$load_per_case <- suggested_load
       },
       required = c(
         ".metainfo$measurements_cols", ".metainfo$flows_cols",
         ".sewer_data$concentrations_observe$measurements",
         ".sewer_data$flows_observe$flows"
         ),
       modeldata = modeldata
       )
     }
     else if (!is.null(min_cases)) {
       modeldata <- tbc("load_per_case_min_cases_calibrate", {
         min_load <- min(modeldata$.metainfo$load_curve_crude$load, na.rm = TRUE)
         load_per_case <- min_load/min_cases
         load_per_case <- signif(load_per_case, signif_fig)
         modeldata$load_mean <- load_per_case
         modeldata$.metainfo$load_per_case <- load_per_case
       },
       required = c(".metainfo$load_curve_crude"),
       modeldata = modeldata
       )
     }
     return(modeldata)
   },
   required_assumptions = "min_cases|cases",
   required_data = "min_cases|cases",
   modeldata = modeldata
  )

  modeldata$.str$shedding[["load_per_case"]] <- list(
    load_per_case_calibrate = c()
  )

  modeldata$.checks$check_load_per_case_function <- function(md, d, a) {
    if ("load_per_case" %in% names(a) && !is.null(a$load_per_case)) {
      cli::cli_inform(c("i"=paste0(
        "You supplied a `load_per_case` assumption via ",
        cli_help("sewer_assumptions"), ", but EpiSewer currently uses ",
        cli_help("load_per_case_calibrate"), " to determine the load per case.",
        " If you want to use the supplied `load_per_case` assumption, please ",
        "specify ", cli_help("load_per_case_assume"), " instead."
      )))
    }
  }

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
  modeldata$nu_zeta_prior <- numeric(0)
  modeldata$.init$nu_zeta <- numeric(0)

  modeldata$zeta_normal_approx <- numeric(0)
  modeldata$n_zeta_normal_approx <- 0
  modeldata$zeta_exact <- numeric(0)
  modeldata$n_zeta_exact <- 0

  modeldata$.init$zeta_log_exact <- numeric(0)
  modeldata$.init$zeta_raw_approx <- numeric(0)

  modeldata$.str$shedding[["load_variation"]] <- list(
    load_variation_none = c()
  )

  return(modeldata)
}

#' Estimate individual-level load variation
#'
#' @description This option accounts for variation in the total shedding load
#'   per case by modeling individual shedding loads as Gamma distributed with
#'   mean equal to the average `load_per_case` and a coefficient of variation to
#'   be estimated.
#'
#' @param cv_prior_mu Mean of the truncated normal prior for the coefficient of
#'   individual-level variation. Default is a CV of 1, which means that
#'   approximately 60% of individuals shed less than the mean, and approximately
#'   10% of individuals shed less than 10% of the mean.
#' @param cv_prior_sigma Standard deviation of the truncated normal prior for
#'   the coefficient of individual-level variation. If this is set to zero, the
#'   coefficient of variation is fixed to the mean.
#'
#' @details Note that measurement noise and individual-level load variation
#'   might not be jointly identifiable. This is why the coefficient of variation
#'   of the individual-level load is fixed by default.
#'
#' @details Also note that the accuracy of the variation model depends on the
#'   estimated number of cases to be on the right scale (which in turn depends
#'   on the assumed `load_per_case` to be roughly correct). This is because in
#'   the variation model, the population-level coefficient of variation (CV) of
#'   shedding loads is proportional to the individual-level CV divided by the
#'   square root of the number of cases. If for example the assumed
#'   `load_per_case` is too small, the number of cases will be overestimated,
#'   which has effects in two directions:
#'   - the individual-level CV will be overestimated (especially if the prior
#'   on the individual-level CV is weak)
#'   - the population-level CV will be underestimated (especially if the prior
#'   on the individual-level CV is strong)
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {load variation models}
load_variation_estimate <- function(
    cv_prior_mu = 1,
    cv_prior_sigma = 0,
    modeldata = modeldata_init()) {
  modeldata$load_vari <- 1
  modeldata$nu_zeta_prior <- set_prior(
    "nu_zeta", "truncated normal",
    mu = cv_prior_mu, sigma = cv_prior_sigma
  )
  modeldata$.init$nu_zeta <- as.array(init_from_location_scale_prior(
    modeldata$nu_zeta_prior
  ))

  modeldata <- tbc("zeta_normal_approximation",
    {
      icc <- modeldata$.metainfo$infection_curve_crude
      l_I <- modeldata$.metainfo$length_I
      l_shedding <- modeldata$.metainfo$length_shedding
      if (nrow(icc) != l_I) {
        cli::cli_abort(
          "Calibration error: Incorrect length of crude infection curve."
        )
      } else {
        infs <- icc[(nrow(icc) - l_shedding + 1):nrow(icc), "infections"][[1]]
        upper_cv <- extraDistr::qtnorm(
          p = 0.95,
          mean = modeldata$nu_zeta_prior$nu_zeta_prior[1],
          sd = modeldata$nu_zeta_prior$nu_zeta_prior[2] + 0.001,
          a = 0
          )
        N_threshold <- 30 * upper_cv^2 # threshold at which to use normal approx
        modeldata$zeta_normal_approx <- which(infs > N_threshold)
        modeldata$n_zeta_normal_approx <- length(modeldata$zeta_normal_approx)
        modeldata$zeta_exact <- which(infs <= N_threshold)
        modeldata$n_zeta_exact <- length(modeldata$zeta_exact)

        modeldata$.init$zeta_log_exact <- rep(
            log(max(1.1, modeldata$.metainfo$initial_cases_crude)),
            modeldata$n_zeta_exact
          )
        modeldata$.init$zeta_raw_approx <- rep(
          0,
          modeldata$n_zeta_normal_approx
        )
      }
    },
    required = c(
      ".metainfo$initial_cases_crude",
      ".metainfo$infection_curve_crude",
      ".metainfo$length_I",
      ".metainfo$length_shedding"),
    modeldata = modeldata
  )

  modeldata$.str$shedding[["load_variation"]] <- list(
    load_variation_estimate = c()
  )

  return(modeldata)
}
