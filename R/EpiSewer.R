#' Estimate the effective reproduction number from wastewater measurements
#'
#' @description `EpiSewer` estimates the effective reproduction number Rt and
#'   other parameters of interest from wastewater measurements over time. This
#'   function combines data, assumptions and modeling details to fit a Bayesian
#'   hierarchical model to the measurements. The resulting posterior estimates
#'   of Rt and other parameters can be plotted and further analyzed.
#'
#' @description The inputs to this function can be specified using various
#' helper functions (see below). It is best to define the inputs in advance as
#' separate variables and then pass them to `EpiSewer`.
#'
#' @param data Observations such as concentration measurements
#' and flows, specified via [sewer_data()]. The data can also be directly
#' supplied to the relevant model components, in which case `data` can be left
#' empty.
#' @param assumptions Assumptions about infections, shedding,
#' or the sewage system, specified via [sewer_assumptions()]. The assumptions
#' can also be directly supplied to the relevant model components, in which
#' case `assumptions` can be left empty.
#' @param measurements The `measurements` module, see [model_measurements()].
#' @param sampling The `sampling` module, see [model_sampling()].
#' @param sewage The `sewage` module, see [model_sewage()].
#' @param shedding The `shedding` module, see [model_shedding()].
#' @param infections The `infections` module, see [model_infections()].
#' @param fit_opts Settings for model fitting, see [set_fit_opts()].
#' @param run_fit If `TRUE` (the default), the model is fitted immediately. If
#' `FALSE`, the EpiSewerJob object is returned without fitting the model (it
#' can be fitted later or even on a different machine).
#'
#' @return If the model fitting was successful, a `list` with the following
#' elements is returned:
#' - `job`: the `EpiSewerJob` that was run
#' - `summary`: a summary of parameter estimates of interest
#' - `fitted`: the fitted model (for `set_fit_opts(fitted = TRUE)`)
#'
#' If the model fitting fails, a `list` with the following elements is returned:
#' - `error`: information about the errors and warnings that were thrown
#' - `sampler_output`: potential outputs printed by the sampler
#'
#' @export
#' @import data.table
#' @import rlang
EpiSewer <- function(
    data = sewer_data(),
    assumptions = sewer_assumptions(),
    measurements = model_measurements(),
    sampling = model_sampling(),
    sewage = model_sewage(),
    shedding = model_shedding(),
    infections = model_infections(),
    fit_opts = set_fit_opts(),
    run_fit = TRUE) {
  modeldata <- modeldata_combine(
    measurements, sampling, sewage, shedding, infections
  )

  model = get_stan_model(modeldata = modeldata)
  modeldata <- modeldata_validate(modeldata, model_def = model)

  job <- EpiSewerJob(
    job_name = paste("Job on", date()),
    model_def = model,
    modeldata = modeldata,
    fit_opts = fit_opts,
    overwrite = TRUE,
    results_exclude = c()
  )

  if (run_fit) {
    return(run(job))
  } else {
    return(list(job = job))
  }
}

#' Constructor for EpiSewerJob objects
EpiSewerJob <- function(job_name,
                        model_def,
                        modeldata,
                        fit_opts,
                        jobarray_size = 1,
                        overwrite = TRUE,
                        results_exclude = c()) {
  job <- list()

  job[["job_name"]] <- job_name
  job[["jobarray_size"]] <- jobarray_size

  job[["model_def"]] <- model_def

  # ToDo rlang::flatten is deprecated, replace
  data_arguments <- suppressWarnings(
    rlang::flatten(modeldata[!(names(modeldata) %in% c("init", "meta_info", "checks"))])
  )
  data_arguments_raw <- data_arguments[
    stringr::str_detect(names(data_arguments), c("_prior_text"), negate = TRUE)
  ]
  job[["data"]] <- data_arguments_raw
  job[["init"]] <- modeldata$init
  job[["fit_opts"]] <- fit_opts

  job[["priors_text"]] <- data_arguments[
    stringr::str_detect(names(data_arguments), "_prior_text")
  ]
  job[["meta_info"]] <- modeldata$meta_info

  job[["overwrite"]] <- overwrite
  job[["results_exclude"]] <- results_exclude

  class(job) <- "EpiSewerJob"

  return(job)
}

setGeneric("run", function(x) UseMethod("run", x))

run.EpiSewerJob <- function(job) {
  arguments <- c(
    list(data = job$data),
    init = function() job$init,
    job$fit_opts$sampler
  )

  result <- list()
  result$job <- job

  fitting_successful <- FALSE
  fit_res <- tryCatch(
    {
      fit_res <- withWarnings(suppress_messages(
        do.call(job$model_def$get_stan_model[[1]]()$sample, arguments),
        "Registered S3 method overwritten by 'data.table'"
      ))
      if (length(fit_res$warnings) == 0) {
        fitting_successful <- TRUE
        fit_res <- fit_res$value
      } else {
        cat("\n")
        rlang::warn(
          paste("There was an error while fitting the model.",
          "Only the model input is returned."))
        fit_res <- list(
          errors = unlist(lapply(fit_res$warnings, function(x) stringr::str_remove(x$message, "\n"))),
          sampler_output = fit_res$value$output()
        )
      }
      fit_res
    },
    error = function(err) {
      cat("\n")
      rlang::warn(c(
        paste("There was an error while fitting the model.",
        "Only the model input is returned."),
        err$message
      ))
      return(list(errors = err, sampler_output = NULL))
    }
  )

  if (!fitting_successful) {
    result$errors <- fit_res$errors
    result$sampler_output <- fit_res$sampler_output
  } else {
    result$summary <- summarize_fit(fit_res, job$data, job$meta_info)

    if (job$fit_opts$fitted) {
      fit_res$draws()
      try(fit_res$sampler_diagnostics(), silent = TRUE)
      try(fit_res$init(), silent = TRUE)
      try(fit_res$profiles(), silent = TRUE)
      result$fitted <- fit_res
    }
  }

  return(result)
}

#' Specify observation data
#'
#' @description Specify wastewater observation data such as concentration
#'   measurements and flows. This is a convenience function to collect all
#'   observation data in one object.
#'
#' @param measurements A `data.frame` with measured concentrations of the
#'   pathogen of interest. Will be automatically passed to
#'   [concentrations_observe()].
#' @param flows A `data.frame` with wastewater flow volumes at the sampling site
#'   for each day. Will be automatically passed to [flows_observe()].
#' @param ... Further observations to be supplied to [EpiSewer()].
#'
#' @return A `list` with all observations supplied. Can be passed to the `data` argument
#'   in [EpiSewer()].
#' @export
sewer_data <- function(measurements = NULL, flows = NULL, ...) {
  data <- as.list(environment())
  return(data)
}

#' Specify modeling assumptions
#'
#' @description Specify model assumptions such as the generation time
#'   distribution or shedding load distribution. This is a convenience function
#'   to collect all assumptions in one object.
#'
#' @param generation_dist Generation time distribution. The intrinsic
#'   distribution of the time between infection of a primary case and infection
#'   of its secondary cases. Will be automatically passed to
#'   [generation_dist_assume()].
#' @param incubation_dist Incubation period distribution. `EpiSewer` uses this
#'   as a proxy for the time between infection and the start of shedding, as
#'   shedding load distributions in the literature are often given from symptom
#'   onset onwards. If the supplied shedding load distribution instead starts
#'   with the time of infection, use `incubation_dist=c(1)` (i.e. no lag). Will
#'   be automatically passed to [incubation_dist_assume()].
#' @param shedding_dist Shedding load distribution. Describes how the total load
#'   shed by an individual is distributed over time (and therefore sums to 1).
#'   Will be automatically passed to [shedding_dist_assume()].
#' @param load_per_case Average total load per case. This is a scaling factor
#'   that describes how many pathogen particles are shed by the average infected
#'   individual overall and how much of this is detectable at the sampling site.
#'   It depends both on biological factors as well as on the specific
#'   sewage system. See [suggest_load_per_case()] to help you
#'   make a suitable assumption. Will be automatically passed to [load_per_case_assume()].
#'
#' @return A `list` with all assumptions supplied. Can be passed to the
#'   `assumptions` argument in [EpiSewer()].
#' @export
sewer_assumptions <- function(generation_dist = NULL,
                              incubation_dist = NULL,
                              shedding_dist = NULL,
                              load_per_case = NULL,
                              residence_dist = c(1),
                              ...) {
  assumptions <- as.list(environment())
  return(assumptions)
}

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
#' `r component_helpers_("concentrations")`
#' @param noise Measurement noise due to unexplained variation in sampling
#' and lab analysis. Modeling options:
#' `r component_helpers_("noise")`
#' @param LOD Limit of detection. Concentrations below a certain threshold may
#' not be detectable and thus erroneously measured as 0. `EpiSewer` can adjust
#' for the limit of detection using a hurdle model. Modeling options:
#' `r component_helpers_("LOD")`
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

#' Model the sampling process
#'
#' @description This module function is used to specify the components of the
#'  `sampling` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param sample_effects Sample (batch) effects. The pathogen concentration in
#' a sample may be influenced by sampling-related external factors, for example
#' the time between sampling and shipping to the lab (age-of-sample effect),
#' or different sampling or storage methods. `EpiSewer` allows to estimate
#' such effects using covariates that describe differences between the samples.
#' Modeling options:
#' `r component_helpers_("sample_effects")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `sampling` module.
#' @export
model_sampling <- function(
    sample_effects = sample_effects_none()) {
  verify_is_modeldata(sample_effects, "sample_effects")
  return(modeldata_combine(sample_effects))
}

#' Model the sewage process
#'
#' @description This module function is used to specify the components of the
#'  `sewage` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param flows Daily flow volumes at the sampling site. The flow can change due
#'   to rainfall or industrial discharge, and directly influences pathogen
#'   concentrations in the wastewater. Modeling options:
#' `r component_helpers_("flows")`
#' @param residence_dist Sewer residence time distribution for pathogen
#'   particles. By default, `EpiSewer` assumes that particles arrive at the
#'   sampling site within the day of shedding. However, for larger sewage
#'   systems, particles may travel longer than a day depending on where and
#'   when they were shed into the wastewater. Modeling options:
#' `r component_helpers_("residence_dist")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `sewage` module.
#' @export
model_sewage <- function(
    flows = flows_observe(),
    residence_dist = residence_dist_assume()) {
  verify_is_modeldata(flows, "flows")
  verify_is_modeldata(residence_dist, "residence_dist")
  return(modeldata_combine(flows, residence_dist))
}

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
#' `r component_helpers_("incubation_dist")`
#' @param shedding_dist Shedding load distribution. Describes how the total load
#'   shed by an individual is distributed over time (and therefore sums to 1).
#'   Modeling options:
#' `r component_helpers_("shedding_dist")`
#' @param load_per_case Average total load per case. This is a scaling factor
#'   that describes how many pathogen particles are shed by the average infected
#'   individual overall and how much of this is detectable at the sampling site.
#'   It depends both on biological factors as well as on the specific
#'   sewage system. Modeling options:
#' `r component_helpers_("load_per_case")`
#' @param load_variation Individual-level shedding load variation. The strength
#' of shedding may vary between individuals. Modeling this variation can better
#' account for uncertainty especially at low incidence. Modeling options:
#' `r component_helpers_("load_variation")`
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
  return(modeldata_combine(incubation_dist, shedding_dist, load_per_case, load_variation))
}

#' Model the infection process
#'
#' @description This module function is used to specify the components of the
#'  `infections` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param generation_dist Generation time distribution. The intrinsic
#'   distribution of the time between infection of a primary case and infection
#'   of its secondary cases. Modeling options:
#' `r component_helpers_("generation_dist")`
#' @param R Effective reproduction number over time. This is the main parameter
#'   of interest estimated by `EpiSewer`. `R` is smoothed using a time series
#'   smoothing prior. Currently supported are: random walk (rw), exponential
#'   smoothing (ets), and smoothing splines. Modeling options:
#' `r component_helpers_("R")`
#' @param seeding Seeding of initial infections. The renewal model used by
#'   `EpiSewer` requires a seeding phase of the length of the maximum generation
#'   time. For these initial infections, a simple seeding model instead of the
#'   renewal model must be used. Modeling options:
#' `r component_helpers_("seeding")`
#' @param infection_noise Noise in the infection process. `EpiSewer` implements
#' a stochastic infection model, i.e. allows for variation in the number of new
#' infections generated at each time step. This accounts for stochastic
#' uncertainty in the infection process and often speeds up model fitting.
#' Modeling options:
#' `r component_helpers_("infection_noise")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `infections` module.
model_infections <- function(
    generation_dist = generation_dist_assume(),
    R = R_estimate_rw(),
    seeding = seeding_estimate(),
    infection_noise = infection_noise_estimate()) {
  verify_is_modeldata(generation_dist, "generation_dist")
  verify_is_modeldata(R, "R")
  verify_is_modeldata(seeding, "seeding")
  verify_is_modeldata(infection_noise, "infection_noise")
  return(modeldata_combine(generation_dist, R, seeding, infection_noise))
}
