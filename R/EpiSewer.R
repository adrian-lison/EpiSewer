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
#'
#' @return
#' @export
#'
#' @examples
EpiSewer <- function(
    data = sewer_data(),
    assumptions = sewer_assumptions(),
    measurements = model_measurements(
      measurements = measurements_observe(measurements = data$measurements)
    ),
    sampling = model_sampling(),
    sewage = model_sewage(
      flows = flows_observe(flows = data$flows)
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
    model = get_stan_model(standata = standata),
    fit_opts = set_fit_opts(sampler_opts = set_sampler_opts(), fitted = TRUE)) {
  standata <- standata_combine(
    measurements, sampling, sewage, shedding, infections
  )
  standata <- standata_validate(standata, model_def = model)


  job <- EpiSewerJob(
    job_name = paste("Job on", date()),
    model_def = model,
    standata = standata,
    fit_opts = fit_opts,
    overwrite = TRUE,
    results_exclude = c()
  )

  fitted <- do.call(job$model_def$get_stan_model[[1]]()$sample, job$arguments)

  result <- list(
    job = job,
    summary = summarize_fit(
      fitted, job$arguments$data, job$arguments_meta_info
    )
  )

  if (fit_opts$fitted) {
    result$fitted = fitted
  }

  return(result)
}

#' Title
#'
#' @param job_name
#' @param model_def
#' @param standata
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
                        standata,
                        fit_opts,
                        jobarray_size = 1,
                        overwrite = TRUE,
                        results_exclude = c()) {
  job_def <- list()

  job_def[["job_name"]] <- job_name
  job_def[["jobarray_size"]] <- jobarray_size

  job_def[["model_def"]] <- model_def

  # ToDo rlang::flatten is deprecated, replace
  data_arguments <- suppressWarnings(
    rlang::flatten(standata[!(names(standata) %in% c("init", "meta_info"))])
  )
  data_arguments_raw <- data_arguments[
    stringr::str_detect(names(data_arguments), c("_prior_text"), negate = TRUE)
  ]
  init_arguments <- standata$init
  job_def[["arguments"]] <- c(list(data = data_arguments_raw),
    init = function() init_arguments, fit_opts$sampler_opts
  )
  job_def[["priors_text"]] <- data_arguments[
    stringr::str_detect(names(data_arguments), "_prior_text")
  ]
  job_def[["arguments_meta_info"]] <- standata$meta_info

  job_def[["overwrite"]] <- overwrite
  job_def[["results_exclude"]] <- results_exclude

  return(job_def)
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
sewer_data <- function(measurements = NULL, flows = NULL) {
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
                              load_per_case = NULL) {
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
    measurements,
    noise = measurement_noise_estimate()) {
  return(standata_combine(measurements, noise))
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
  return(standata_combine(sample_effects))
}

#' Title
#'
#' @param flows
#'
#' @return
#' @export
#'
#' @examples
model_sewage <- function(flows) {
  return(standata_combine(flows))
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
model_shedding <- function(incubation_dist, shedding_dist, load_per_case) {
  return(standata_combine(incubation_dist, shedding_dist, load_per_case))
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
    generation_dist,
    R = R_estimate_splines(),
    seeding = seeding_estimate(),
    infection_noise = infection_noise_estimate()) {
  return(standata_combine(generation_dist, R, seeding, infection_noise))
}
