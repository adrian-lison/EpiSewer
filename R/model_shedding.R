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
model_shedding <- function(
    incubation_dist = incubation_dist_assume(),
    shedding_dist = shedding_dist_assume(),
    load_per_case = load_per_case_assume(),
    load_variation = load_variation_none()) {
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
#' @details `EpiSewer` uses the incubation period as a proxy for the time
#'   between infection and the start of shedding. This is because shedding load
#'   distributions in the literature are often given from symptom onset onwards.
#'   If the assumed shedding load distribution instead starts from the time of
#'   infection, then use `incubation_dist=c(1)` (i.e. no lag).
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
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso Helpers to discretize continuous probability distributions:
#'   [get_discrete_gamma()], [get_discrete_lognormal()]
shedding_dist_assume <-
  function(shedding_dist = NULL, modeldata = modeldata_init()) {
    modeldata <- tbp("shedding_dist_assume",
      {
        modeldata$S <- length(shedding_dist) - 1
        shedding_dist <- check_dist(shedding_dist, "shedding load distribution")
        modeldata$shedding_dist <- shedding_dist
        return(modeldata)
      },
      required_assumptions = "shedding_dist",
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
#' @seealso {Helper for finding a suitable load per case assumption}
#'   [suggest_load_per_case()]
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
  modeldata$.init$zeta <- numeric(0)

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
#' @details Measurement noise and individual-level load variation might not be
#'   jointly identifiable. It may therefore be necessary to use informative
#'   priors for at least one of these parameters when using both
#'   [noise_estimate()] and [load_variation_estimate()].
#'
#' @details Also note that the accuracy of the variation model depends on
#'   the estimated number of cases to be on the right scale (which in turn
#'   depends on the assumed `load_per_case` to be roughly correct). This is
#'   because in the variation model, the population-level coefficient of
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
  modeldata$nu_zeta_prior <- set_prior(
    "nu", "truncated normal",
    mu = cv_prior_mu, sigma = cv_prior_sigma
  )
  modeldata$.init$nu_zeta <- as.array(cv_prior_mu)
  modeldata$.init$zeta <- tbe(
    rep(
      median(modeldata$measured_concentrations, na.rm = T),
      modeldata$.metainfo$length_shedding
    ),
    c("measured_concentrations", ".metainfo$length_shedding")
  )

  modeldata$.str$shedding[["load_variation"]] <- list(
    load_variation_estimate = c()
  )

  return(modeldata)
}
