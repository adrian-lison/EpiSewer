#' Title
#'
#' @param data
#' @param assumptions
#' @param measurements
#' @param sampling
#' @param sewage
#' @param shedding
#' @param infections
#' @param model
#' @param fit_opts
#' @param run_fit Run the model fitting immediately or return only the job without running it?
#'
#' @return
#' @export
#' @import data.table
#' @import rlang
#'
#' @examples
EpiSewer <- function(
    data = sewer_data(),
    assumptions = sewer_assumptions(),
    measurements = model_measurements(
      concentrations = concentrations_observe(data = data$measurements)
    ),
    sampling = model_sampling(),
    sewage = model_sewage(
      flows = flows_observe(data = data$flows),
      residence_dist = residence_dist_assume(
        residence_dist = assumptions$residence_dist
        )
    ),
    shedding = model_shedding(
      incubation_dist = incubation_dist_assume(
        incubation_dist = assumptions$incubation_dist
      ),
      shedding_dist = shedding_dist_assume(
        shedding_dist = assumptions$shedding_dist
      ),
      load_per_case = load_per_case_assume(
        load_per_case = assumptions$load_per_case
      )
    ),
    infections = model_infections(
      generation_dist = generation_dist_assume(
        generation_dist = assumptions$generation_dist
      )
    ),
    model = get_stan_model(modeldata = modeldata),
    fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(), fitted = TRUE),
    run_fit = TRUE) {
  modeldata <- modeldata_combine(
    measurements, sampling, sewage, shedding, infections
  )
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

#' Title
#'
#' @param job_name
#' @param model_def
#' @param modeldata
#' @param fit_opts
#' @param jobarray_size
#' @param overwrite
#' @param results_exclude
#'
#' @return
#' @export
#'
#' @examples
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
        rlang::warn(c(
          "There was an error while fitting the model:",
          "Only the model input is returned."))
        fit_res <- unlist(lapply(fit_res$warnings, function(x) stringr::str_remove(x$message, "\n")))
      }
      fit_res
    },
    error = function(err) {
      cat("\n")
      rlang::warn(c(
        "There was an error while fitting the model:",
        "Only the model input is returned.",
        err$message
      ))
      return(err)
    }
  )

  if (!fitting_successful) {
    result$errors <- fit_res
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

#' Title
#'
#' @param measurements
#' @param flows
#'
#' @return
#' @export
#'
#' @examples
sewer_data <- function(measurements = NULL, flows = NULL, ...) {
  data <- as.list(environment())
  return(data)
}

#' Title
#'
#' @param generation_dist
#' @param incubation_dist
#' @param shedding_dist
#' @param load_per_case
#'
#' @return
#' @export
#'
#' @examples
sewer_assumptions <- function(generation_dist = NULL,
                              incubation_dist = NULL,
                              shedding_dist = NULL,
                              load_per_case = NULL,
                              residence_dist = c(1),
                              ...) {
  assumptions <- as.list(environment())
  return(assumptions)
}

#' Title
#'
#' @param measurements
#' @param noise
#'
#' @return
#' @export
#'
#' @examples
model_measurements <- function(
    concentrations = concentrations_observe(),
    noise = noise_estimate()) {
  verify_is_modeldata(concentrations, "concentrations")
  verify_is_modeldata(noise, "noise")
  return(modeldata_combine(concentrations, noise))
}

#' Title
#'
#' @param sample_effects
#'
#' @return
#' @export
#'
#' @examples
model_sampling <- function(
    sample_effects = sample_effects_none()) {
  verify_is_modeldata(sample_effects, "sample_effects")
  return(modeldata_combine(sample_effects))
}

#' Title
#'
#' @param flows
#'
#' @return
#' @export
#'
#' @examples
model_sewage <- function(
    flows = flows_observe(),
    residence_dist = residence_dist_assume()) {
  verify_is_modeldata(flows, "flows")
  verify_is_modeldata(residence_dist, "residence_dist")
  return(modeldata_combine(flows, residence_dist))
}

#' Title
#'
#' @param incubation_dist
#' @param shedding_dist
#' @param load_per_case
#'
#' @return
#' @export
#'
#' @examples
model_shedding <- function(
    incubation_dist = incubation_dist_assume(),
    shedding_dist = shedding_dist_assume(),
    load_per_case = load_per_case_assume()) {
  verify_is_modeldata(incubation_dist, "incubation_dist")
  verify_is_modeldata(shedding_dist, "shedding_dist")
  verify_is_modeldata(load_per_case, "load_per_case")
  return(modeldata_combine(incubation_dist, shedding_dist, load_per_case))
}

#' Title
#'
#' @param generation_dist
#' @param R
#' @param seeding
#' @param infection_noise
#'
#' @return
#' @export
#'
#' @examples
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
