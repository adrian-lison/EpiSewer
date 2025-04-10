#' Estimate epidemiological parameters from wastewater measurements
#'
#' @description `EpiSewer` estimates the effective reproduction number Rt and
#'   other parameters of interest from wastewater measurements over time. This
#'   function combines data, assumptions and modeling details to fit a Bayesian
#'   hierarchical model to the measurements. The resulting posterior estimates
#'   of Rt and other parameters can be plotted and further analyzed.
#'
#' @description The inputs to this function can be specified using various
#'   helper functions (see below). It is best to define the inputs in advance as
#'   separate variables and then pass them to `EpiSewer`.
#'
#' @param data Observations such as concentration measurements and flows,
#'   specified via [sewer_data()]. The data can also be directly supplied to the
#'   relevant model components, in which case `data` can be left empty.
#' @param assumptions Assumptions about infections, shedding, or the sewage
#'   system, specified via [sewer_assumptions()]. The assumptions can also be
#'   directly supplied to the relevant model components, in which case
#'   `assumptions` can be left empty.
#' @param measurements The `measurements` module, see [model_measurements()].
#' @param sampling The `sampling` module, see [model_sampling()].
#' @param sewage The `sewage` module, see [model_sewage()].
#' @param shedding The `shedding` module, see [model_shedding()].
#' @param infections The `infections` module, see [model_infections()].
#' @param fit_opts Settings for model fitting, see [set_fit_opts()].
#' @param results_opts Settings for results to be returned, see
#'   [set_results_opts()].
#' @param run_fit If `TRUE` (the default), the model is fitted immediately. If
#'   `FALSE`, the EpiSewerJob object is returned without fitting the model (it
#'   can be fitted later or even on a different machine).
#'
#' @return If the model fitting was successful, a `list` with the following
#'   elements is returned:
#' - `job`: the `EpiSewerJob` that was run
#' - `summary`: a summary of parameter estimates of interest
#' - `fitted`: the fitted model (if `results_opts(fitted = TRUE)`)
#'
#' If the model fitting fails, a `list` with the following elements is
#'   returned:
#' - `error`: information about the errors and warnings that were thrown
#' - `sampler_output`: potential outputs printed by the sampler
#'
#' @export
#' @import data.table
EpiSewer <- function(
    data = sewer_data(),
    assumptions = sewer_assumptions(),
    measurements = model_measurements(),
    sampling = model_sampling(),
    sewage = model_sewage(),
    shedding = model_shedding(),
    infections = model_infections(),
    forecast = model_forecast(),
    fit_opts = set_fit_opts(),
    results_opts = set_results_opts(),
    run_fit = TRUE) {
  modeldata <- modeldata_combine(
    measurements, sampling, sewage, shedding, infections, forecast
  )

  modeldata <- modeldata_validate(
    modeldata, data = data, assumptions = assumptions
  )

  job <- EpiSewerJob(
    job_name = paste("EpiSewerJob", format(lubridate::now(), "%Y-%m-%d_%H-%M-%S")),
    modeldata = modeldata,
    fit_opts = fit_opts,
    results_opts = results_opts,
    overwrite = TRUE,
    results_exclude = c()
  )

  if (run_fit) {
    return(run(job))
  } else {
    res <- list(job = job)
    class(res) <- c("EpiSewerJobResult", class(res))
    return(res)
  }
}

#' Constructor for EpiSewerJob objects
#' @keywords internal
EpiSewerJob <- function(job_name,
                        modeldata,
                        fit_opts,
                        results_opts,
                        jobarray_size = 1,
                        overwrite = TRUE,
                        results_exclude = c()) {
  job <- list()

  job[["job_name"]] <- job_name
  job[["jobarray_size"]] <- jobarray_size

  # ToDo rlang::flatten is deprecated, replace
  data_arguments <- suppressWarnings(
    rlang::flatten(modeldata[!(names(modeldata) %in% c(
      ".init", ".metainfo", ".checks", ".str",
      ".sewer_data", ".sewer_assumptions"
    ))])
  )
  data_arguments_raw <- data_arguments[
    stringr::str_detect(names(data_arguments), c("_prior_text"), negate = TRUE)
  ]
  job[["data"]] <- data_arguments_raw
  job[["model"]] <- modeldata$.str
  job[["init"]] <- modeldata$.init
  job[["fit_opts"]] <- fit_opts
  job[["results_opts"]] <- results_opts

  job[["priors_text"]] <- data_arguments[
    stringr::str_detect(names(data_arguments), "_prior_text")
  ]
  job[["metainfo"]] <- modeldata$.metainfo

  job[["overwrite"]] <- overwrite
  job[["results_exclude"]] <- results_exclude

  class(job) <- c("EpiSewerJob", class(job))

  return(job)
}

#' @export
setClass("EpiSewerJob")

#' Fit an EpiSewer model.
#'
#' @description This function allows to (re-)run the model fitting from an
#'   existing EpiSewerJob object. It is especially useful when combined with the
#'   `run_fit=FALSE` option in [EpiSewer()]. This allows to first create the
#'   EpiSewerJob object and then run the model fitting later or on a different
#'   machine.
#'
#' @param job An EpiSewerJob object as returned by [EpiSewer()].
#'
#' @export
run <- function(job) {
  UseMethod("run")
}

#' @export
run.EpiSewerJob <- function(job) {
  arguments <- c(
    list(data = job$data),
    init = function() job$init,
    job$fit_opts$sampler
  )

  fitting_successful <- FALSE
  result <- list()
  result$job <- job

  result$stan_model <- get_stan_model(
    model_metainfo = job$metainfo,
    model_folder = job$fit_opts$model$model_folder,
    profile = job$fit_opts$model$profile,
    threads = job$fit_opts$model$threads,
    force_recompile = job$fit_opts$model$force_recompile,
    package = job$fit_opts$model$package
  )

  stanmodel_instance <- result$stan_model$load_model[[1]]()

  result$checksums <- get_checksums(job, stanmodel_instance)

  fit_res <- tryCatch(
    {
      fit_res <- withWarnings(suppress_messages_warnings(
        do.call(stanmodel_instance$sample, arguments),
        c(
          "Registered S3 method overwritten by 'data.table'",
          "Cannot parse stat file, cannot read file: No such file or directory",
          "cannot open file '/proc/stat': No such file or directory"
        )
      ))
      if (length(fit_res$warnings) == 0) {
        fitting_successful <- TRUE
        fit_res <- fit_res$value
      } else {
        cat("\n")
        cli::cli_warn(
          paste(
            "There was an error while fitting the model.",
            "Only the model input is returned."
          )
        )
        fit_res <- list(
          errors = unlist(lapply(
            fit_res$warnings, function(x) stringr::str_remove(x$message, "\n")
          )),
          sampler_output = fit_res$value$output()
        )
      }
      fit_res
    },
    error = function(err) {
      cat("\n")
      cli::cli_warn(c(
        paste(
          "There was an error while fitting the model.",
          "Only the model input is returned."
        ),
        err$message
      ))
      return(list(errors = err, sampler_output = NULL))
    }
  )

  if (fitting_successful) {
    result$summary <- try(summarize_fit(
      fit = fit_res,
      data = job$data,
      .metainfo = job$metainfo,
      intervals = job$results_opts$summary_intervals,
      ndraws = job$results_opts$samples_ndraws
      ))
    if (job$results_opts$fitted) {
      try(fit_res$draws(), silent = TRUE)
      try(fit_res$sampler_diagnostics(), silent = TRUE)
      try(fit_res$init(), silent = TRUE)
      try(fit_res$profiles(), silent = TRUE)
      result$fitted <- fit_res
    }
    result$diagnostics <- try(suppressMessages(fit_res$diagnostic_summary()))
    result$runtime <- try(fit_res$time())
  } else {
    result$errors <- fit_res$errors
    result$sampler_output <- fit_res$sampler_output
  }

  class(result) <- c("EpiSewerJobResult", class(result))
  return(result)
}

#' @export
run.EpiSewerJobResult <- function(job) {
  return(run(job$job))
}

#' @export
test_run <- function(job) {
  UseMethod("test_run")
}

#' @export
test_run.EpiSewerJob <- function(job) {
  job$fit_opts$sampler$iter_warmup <- 5
  job$fit_opts$sampler$iter_sampling <- 5
  job$fit_opts$sampler$show_messages <- FALSE
  job$fit_opts$sampler$show_exceptions <- FALSE
  job$fit_opts$sampler$output_dir <- withr::local_tempdir()
  job$results_opts$samples_ndraws <- 2
  return(run(job))
}

#' @export
test_run.EpiSewerJobResult <- function(job) {
  return(test_run(job$job))
}

#' Print an EpiSewerJob.
#' @export
#' @keywords internal
print.EpiSewerJob <- function(x) {
  cat(paste0("\n",x$job_name,"\n"))
  cat("-----\n")
  cat(paste0(
    "Time period: ",
    x$metainfo$T_start_date, " to ", x$metainfo$T_end_date, "\n",
    "Observations: ", x$data$n_measured, " measurements on ",
    length(unique(x$data$sample_to_date)), " days\n"
    ))
  if (x$metainfo$forecast_horizon > 0) {
    cat(paste0(
      "Forecast horizon: ", x$metainfo$forecast_horizon, " days\n"
    ))
  }
  cat("\nUse ...$model for a summary of the model specification.")
}

#' @export
setClass("EpiSewerJobResult")

#' Print an EpiSewerJobResult.
#' @export
#' @keywords internal
print.EpiSewerJobResult <- function(x) {
  print(x$job)
  if (all(c("runtime","diagnostics") %in% names(x))) {
    cat("\n\n")
    cat(paste0("Model fitted in ", round(x$runtime$total),
               " seconds (", nrow(x$runtime$chains) ," chains, ",
               x$job$fit_opts$sampler$iter_sampling, " iterations each)\n"))
    cat("-----")
    num_divergent <- sum(x$diagnostics$num_divergent)
    if (num_divergent>0) cat("\nDivergent transitions:", num_divergent)
    num_max_treedepth <- sum(x$diagnostics$num_max_treedepth)
    if (num_max_treedepth>0) cat("\nMaximum treedepth reached:", num_max_treedepth)
    num_ebfmi_low <- sum(x$diagnostics$num_ebfmi_low<0.2)
    if (num_ebfmi_low>0) cat("\nE-BFMI problems:", num_ebfmi_low)
    if (num_divergent==0 && num_max_treedepth==0 && num_ebfmi_low==0) {
      cat("\nAll diagnostics are perfect.")
    }
  } else if ("errors" %in% names(x)) {
    cat("\n\n")
    cat("Model fitting failed:\n  ")
    cat(paste(x$errors, collapse="\n  "))
    cat("\nDetailed sampler output can be accessed via ...$sampler_output.")
  }
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
#' @param cases A `data.frame` of case numbers with each row representing one
#'   day. Must have at least a column with dates and a column with case numbers.
#'   This data is not used for model fitting, but for calibration of the
#'   `load_per_case` assumption. Will be automatically passed to
#'   [load_per_case_calibrate()].
#' @param ... Further observations to be passed to [EpiSewer()].
#'
#' @return A `list` with all observations supplied. Can be passed to the `data`
#'   argument in [EpiSewer()].
#' @export
sewer_data <- function(measurements = NULL, flows = NULL, cases = NULL, ...) {
  data <- c(as.list(environment()), list(...))
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
#' @param shedding_dist Shedding load distribution. Describes how the total load
#'   shed by an individual is distributed over time (and therefore sums to 1).
#'   Will be automatically passed to [shedding_dist_assume()].
#' @param shedding_reference Is the shedding load distribution relative to the
#'   day of `"infection"` or the day of `"symptom_onset"`? This is important because
#'   shedding load distributions provided in the literature are sometimes by
#'   days since infection and sometimes by days since symptom onset. If
#'   `shedding_reference="symptom_onset"`, EpiSewer also needs information about
#'   the incubation period distribution. Will be automatically passed to
#'   [shedding_dist_assume()].
#' @param incubation_dist Incubation period distribution. The incubation period
#'   is the time between infection and symptom onset. This assumption is used
#'   when `shedding_reference="symptom_onset"`, i.e. to support shedding load
#'   distributions referenced by days since symptom onset. Will be automatically
#'   passed to [incubation_dist_assume()].
#' @param min_cases The minimum incidence (new cases per day) during the modeled
#'   time period (default is 10). This assumption is used to calibrate the
#'   `load_per_case` factor (a scaling factor that describes how many pathogen
#'   particles are shed by the average infected individual overall and how much
#'   of this is detectable at the sampling site). If `min_cases` is specified,
#'   `load_per_case` is chosen such that the estimated time series of infections
#'   does not go significantly below `min_cases`. Will be automatically passed
#'   to [load_per_case_calibrate()].
#' @param load_per_case This argument allows to directly specify the average
#'   total load per case, instead of calibrating it to case data. Note that the
#'   `load_per_case` depends both on biological factors as well as on the
#'   specific sewage system and laboratory quantification. Important: To
#'   actually use this assumption, the model option [load_per_case_assume()]
#'   must be specified in [EpiSewer()] via [model_shedding()]. If you want to
#'   use this option to do a manual calibration to case data, the helper
#'   function [suggest_load_per_case()] might also be useful.
#' @param residence_dist Sewer residence time distribution for pathogen
#'   particles. By default, `EpiSewer` assumes that particles arrive at the
#'   sampling site within the day of shedding. However, for larger sewage
#'   systems, particles may travel longer than a day depending on where and when
#'   they were shed into the wastewater. Will be automatically passed to
#'   [residence_dist_assume()].
#' @param ... Further assumptions to be passed to [EpiSewer()].

#' @details Note that if case data is supplied via `sewer_data`, the `min_cases`
#'   assumption is overwritten and the supplied case data is used for a more
#'   accurate calibration instead.
#'
#' @return A `list` with all assumptions supplied. Can be passed to the
#'   `assumptions` argument in [EpiSewer()].
#' @export
sewer_assumptions <- function(generation_dist = NULL,
                              shedding_dist = NULL,
                              shedding_reference = NULL,
                              incubation_dist = NULL,
                              min_cases = 10,
                              load_per_case = NULL,
                              residence_dist = c(1),
                              ...) {
  assumptions <- c(as.list(environment()), list(...))
  return(assumptions)
}

#' Configure results returned after model fitting
#'
#' @param fitted If `TRUE` (default), the fitted model object is also returned,
#'   not only summaries of the model fitting. Note that this makes the results
#'   object substantially larger.
#' @param summary_intervals Which credible intervals (CrIs) should be used to
#'   summarize the posterior distributions? Default is `c(0.5, 0.95)`, i.e. the
#'   50% and 95% CrI are returned
#' @param samples_ndraws Number of exemplary posterior samples to return. Note
#'   that the summaries always use all available samples.
#'
#' @return A `list` with settings for the results.
#' @export
#' @examples
#' ww_data <- ww_data_SARS_CoV_2_Zurich
#' ww_assumptions <- ww_assumptions_SARS_CoV_2_Zurich
#'
#' ww_fit_opts <- set_fit_opts(
#'   sampler = sampler_stan_mcmc(
#'     chains = 2,
#'     parallel_chains = 2,
#'     iter_warmup = 400,
#'     iter_sampling = 400,
#'     seed = 42 # ensures reproducibility
#'   )
#' )
#'
#' ww_results_opts <- set_results_opts(
#'   fitted = TRUE, # return full model fit
#'   summary_intervals = c(0.5, 0.8, 0.95),
#'   samples_ndraws = 100
#' )
#'
#' \dontrun{ww_result <- EpiSewer(
#'   data = ww_data,
#'   assumptions = ww_assumptions,
#'   fit_opts = ww_fit_opts,
#'   results_opts = ww_results_opts
#' )}
set_results_opts <- function(fitted = TRUE, summary_intervals = c(0.5, 0.95), samples_ndraws = 50) {
  opts <- as.list(environment())
  return(opts)
}

#' Get checksums that uniquely identify an EpiSewer job.
#'
#' @param job An EpiSewer job object.
#' @param stanmodel_instance A stanmodel instance.
#'
#' @return A `list` with checksums identifying the model, input, inits,
#'   fit_opts and results_opts.
#' @export
get_checksums <- function(job, stanmodel_instance = NULL) {
  checksums <- list()
  if (!is.null(stanmodel_instance)) {
    checksums$model <- get_checksum_model(stanmodel_instance)
  } else {
    checksums$model <- NA
  }
  checksums$input <- digest::digest(job$data, algo = "md5")
  checksums$fit_opts <- digest::digest(
    list_except(job$fit_opts, c("model")), algo = "md5"
    )
  checksums$results_opts <- digest::digest(
    job$results_opts, algo = "md5"
  )
  checksums$init <- digest::digest(job$init, algo = "md5")
  return(checksums)
}
