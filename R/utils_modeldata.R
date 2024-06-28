#' Template function
#'
#' @description This function does not actually do anything. It only serves as a
#'   template for the documentation of other functions using inheritance of
#'   params.
#'
#' @param modeldata A `modeldata` object to which the above model specifications
#'   should be added. Default is an empty model given by [modeldata_init()]. Can
#'   also be an already partly specified model returned by other `EpiSewer`
#'   modeling functions.
#'
#' @return Nothing
template_model_helpers <- function(modeldata) { }

#' Construct an unspecified EpiSewer model
#'
#' @return A `modeldata` object containing data and specifications of the model
#'   to be fitted. Can be passed on to other `EpiSewer` modeling functions to
#'   add further data and model specifications.
#'
#' @return The `modeldata` object also includes information about parameter
#'   initialization (`init`), meta data (`.metainfo`), and checks to be
#'   performed before model fitting (`.checks`).
#' @export
#'
#' @examples
#' modeldata_init()
modeldata_init <- function() {
  modeldata <- list()
  modeldata$.init <- list()
  modeldata$.metainfo <- list()
  modeldata$.checks <- list()
  modeldata$.str <- list(
    measurements = list(),
    sampling = list(),
    sewage = list(),
    shedding = list(),
    infections = list()
  )
  class(modeldata$.str) <- c("modelstructure")
  modeldata$.sewer_data <- list()
  modeldata$.sewer_assumptions <- list()
  class(modeldata) <- "modeldata"
  return(modeldata)
}

#' Verification of modeldata objects
verify_is_modeldata <- function(modeldata, arg_name) {
  if (!all(c(".init", ".metainfo", ".checks") %in% names(modeldata))) {
    error_msg <- paste0(
      "The argument `", arg_name, "` expects EpiSewer model data as input. ",
      "Please use a function of the form `", arg_name, "_`", " to specify ",
      "this argument."
    )
    all_fs <- names(rlang::ns_env("EpiSewer"))
    arg_fs <- all_fs[stringr::str_detect(
      # this excludes functions ending with _
      all_fs, paste0("^", arg_name, "_", ".*(?<!_)$")
    )]
    if (length(arg_fs) > 0) {
      functions <- paste0(arg_fs, "()")
      functions_cli <- paste0(
        "{.help [", functions, "](EpiSewer::", functions, "}"
      )
      error_msg <- paste(error_msg, "Available functions:")
      error_msg <- c(error_msg, functions_cli)
      names(error_msg) <- c("!", rep("*", length(error_msg) - 1))
    }
    cli::cli_abort(error_msg, call = rlang::caller_env())
  }
}

#' Get modeling functions for a component (internal function)
#'
#' Get all available modeling functions for a given component. Argument defaults
#' are such that helpers are returned as a comma-separated list of roxygen
#' function links.
#'
#' @param component A character vector with the name of the component
#' @param collapse A character used to collapse helpers into one string. Default
#'   is NULL (do not collapse).
#' @param prefix A prefix to add to each function name
#' @param suffix A suffix to add to each function name
#'
#' @return A character vector with all available modeling functions for the
#'   component. If collapse is not NULL, the helper functions are collapsed into
#'   a single string.
component_functions_ <- function(component, collapse = "\n",
                                 prefix = "- [", suffix = "()]") {
  if (!component %in% all_components()) {
    err_message <- c(
      paste(
        "No valid component provided.",
        "Must be one out of:"
      ),
      all_components()
    )
    names(err_message) <- c("!", rep("*", length(err_message) - 1))
    cli::cli_abort(err_message)
  }
  all_fs <- names(rlang::ns_env("EpiSewer"))
  arg_fs <- all_fs[stringr::str_detect(
    # this excludes functions ending with _
    all_fs, paste0("^", component, "_", ".*(?<!_)$")
  )]

  # sort values in arg_fs according to following scheme
  # 1. values ending with "_none" first
  # 2. values ending with "_observe" second
  # 3. values ending with "_assume" third
  # 4. values ending with "_estimate" fourth
  # 5. values ending with "_estimate" plus another suffix fifth
  arg_fs_none <- arg_fs[stringr::str_detect(arg_fs, "_none$")]
  arg_fs_observe <- arg_fs[stringr::str_detect(arg_fs, "_observe$")]
  arg_fs_assume <- arg_fs[stringr::str_detect(arg_fs, "_assume$")]
  arg_fs_estimate <- arg_fs[stringr::str_detect(arg_fs, "_estimate$")]
  arg_fs_estimate_other <- arg_fs[stringr::str_detect(arg_fs, "_estimate_")]
  arg_fs_estimate_other <- arg_fs_estimate_other[
    !arg_fs_estimate_other %in% c(
      arg_fs_none, arg_fs_observe, arg_fs_assume, arg_fs_estimate
    )
  ]
  arg_fs_sorted <- c(
    arg_fs_none, arg_fs_observe, arg_fs_assume, arg_fs_estimate, arg_fs_estimate_other
  )

  if (length(arg_fs_sorted) > 0) {
    helpers <- paste0(prefix, arg_fs_sorted, suffix)
  } else {
    helpers <- c()
  }
  if (!is.null(collapse)) {
    return(paste(helpers, collapse = collapse))
  } else {
    return(helpers)
  }
}

#' Get modeling functions for a component
#'
#' @inheritParams component_functions_
#'
#' @return A character vector with all available modeling functions for the
#'   component.
#' @export
#'
#' @examples
#' component_functions("R")
component_functions <- function(component) {
  return(component_functions_(
    component,
    collapse = NULL, prefix = "", suffix = "()"
  ))
}

#' Define a prior in modeldata
set_prior <- function(param, dist = "normal", ...) {
  prior <- c(as.list(environment()), list(...))
  prior_data <- list()
  prior_data[[paste0(param, "_prior_text")]] <- paste(
    names(prior), "=", prior,
    collapse = ", "
  )
  prior$param <- NULL
  prior$dist <- NULL
  names(prior) <- NULL
  prior_data[[paste0(param, "_prior")]] <- unlist(prior)
  return(prior_data)
}


#' Define a normal prior in modeldata
#'
#' @param param Name of the parameter for which the prior is defined.
#' @param mu Mean.
#' @param sigma Standard deviation.
#' @param two_sigma Two times the standard deviation. Useful for defining priors
#'   via the two-sigma rule-of-thumb (approximately 95% of probability mass).
#'
#' @return Prior specification for modeldata.
set_prior_normal <- function(param, mu, sigma = NULL, two_sigma = NULL) {
  if (is.null(sigma)) {
    if (is.null(two_sigma)) {
      cli::cli_abort(
        "Either sigma or two_sigma must be supplied for the normal prior",
        .internal = TRUE
        )
    } else {
      sigma <- two_sigma / 2
    }
  }
  return(set_prior(param = param, dist = "normal", mu = mu, sigma = sigma))
}

#' Define a normal prior on the log scale in modeldata using median and factor
#'
#' @description This parameterization is useful to define a best guess for the
#'   value of a parameter (median), and a maximum factor by which we expect to
#'   deviate from that guess.
#'
#' @param param Name of the parameter for which the prior is defined.
#' @param unit_median Median on the unit/natural scale.
#' @param unit_q5 5% quantile of the distribution on the unit scale.
#' @param unit_q95 95% quantile of the distribution on the unit scale.
#' @param unit_factor By which factor do we expect the true parameter value to
#'   differ at most from our prior median? For example, `unit_factor = 2` would
#'   mean that we expect the true parameter to be at most twice our prior
#'   median, and at least half of our prior median. This uses the two-sigma
#'   rule-of-thumb. Must be specified together with `unit_median.`
#' @param paste_dist Additional information that should be pasted to the dist
#'   description.
#'
#' @details A normal prior on the log scale is effectively a log-normal prior on
#'   the unit/natural scale.
#'
#' @return Prior specification for modeldata.
set_prior_normal_log <- function(param,
                                 unit_median = NULL,
                                 unit_q5 = NULL, unit_q95 = NULL,
                                 unit_factor = NULL, paste_dist = "") {
  if (!is.null(unit_median)) {
    mu = log(unit_median)
    if (!is.null(unit_q5) || !is.null(unit_q95)) {
      sigma = get_lognormal_sigma_alternative(
        mu = mu, unit_q5 = unit_q5, unit_q95 = unit_q95
        )
    } else if (!is.null(unit_factor)) {
      sigma = log(unit_factor)/2
    } else {
      cli::cli_abort(paste(
        "You need to specify one additional argument besides",
        "`unit_median` to define the prior."
        ))
    }
  } else {
    if (!is.null(unit_q5) && !is.null(unit_q95)) {
      mu = get_lognormal_mu_alternative(unit_q5 = unit_q5, unit_q95 = unit_q95)
      sigma = get_lognormal_sigma_alternative(mu = mu, unit_q95 = unit_q95)
    } else {
      cli::cli_abort(paste(
        "Please specify either `unit_median`, or both `unit_q5` and",
        "`unit_q95` to define the prior."
      ))
    }
  }
  prior = set_prior(
    param = param, dist = paste0("normal", paste_dist), mu = mu , sigma = sigma
    )
  return(prior)
}

#' Provide initialization value for a parameter based on the supplied prior
#'
#' @description Initialization using the prior is often better than initializing
#'   with zero (and if the parameter is strictly positive, zero is not possible
#'   at all). This function provides as init value the mean of the normal prior
#'   plus 1/4 of the standard deviation. This ensure a positive init even if the
#'   mean is zero (useful for truncated normal priors.)
#'
#' @param prior Prior for parameter as provided by [set_prior()]. Should be a
#'   normal prior (first element mean, second element sd).
#'
#' @details If the provided prior has zero variance, it is assumed that the
#'   parameter will not be sampled and an empty init is returned.
#'
#' @return Mean of the normal prior plus 1/4 of the standard deviation.
init_from_normal_prior <- function(prior) {
  prior_data_select <- !stringr::str_detect(names(prior), pattern = "_text$")
  if (sum(prior_data_select) != 1) {
    cli::cli_warn(paste(
      "Could not init from normal prior:",
      "Non-ambiguous prior format. Using 1e-2 as fallback."
      ), .internal = TRUE)
  }
  prior_data <- prior[prior_data_select][[1]]
  if (prior_data[2] > 0) {
    return(prior_data[1] + prior_data[2]/4)
  } else {
    return(numeric(0))
  }
}

#' Check modeldata object for required variables
#'
#' `modeldata_check` accepts various additional arguments which are turned into
#' a sophisticated error message about missing variables
#'
#' @inheritParams modeldata_validate
#' @param required Character vector of required variables
#' @param descriptions Optional, a list with names corresponding to all or some
#'   of the required variables, and values corresponding to additional
#'   descriptions of the variables
#' @param run_before Optional, either a character vector, or a `list` with names
#'   corresponding to all of the required variables, and values corresponding to
#'   functions that will add these variables to modeldata.
#' @param advice Additional message with advice on how to solve error.
#' @param throw_error Should an error be thrown (with verbatim error message),
#'   or should `TRUE`/`FALSE` be returned. Default is `TRUE`.
#'
#' @return If throw_error is `TRUE`, nothing is returned, otherwise a logical is
#'   returned (`TRUE` if all required variables are present)
modeldata_check <- function(modeldata,
                            required,
                            required_values = NULL,
                            descriptions = modeldata_descriptions(),
                            run_before = modeldata_var_requirements(),
                            advice = NULL,
                            throw_error = TRUE,
                            calling_env = rlang::caller_env()) {
  if (!class(modeldata) == "modeldata") {
    cli::cli_abort("Please supply a modeldata object.")
  }

  var_check <- check_list_nested(modeldata, required, required_values)
  if (throw_error) {
    if (!all(var_check)) {
      missing_vars <- required[!var_check]
      # add descriptions of vars in brackets
      missing_vars_text <- sapply(missing_vars, function(x) {
        if (x %in% names(descriptions)) {
          return(paste0(x, " (", descriptions[[x]], ")"))
        } else {
          return(x)
        }
      })
      error_msg <- paste0(
        "The following variables are required ",
        "but currently not present in modeldata:\n\n",
        paste(paste(" -", missing_vars_text), collapse = "\n"), "\n"
      )
      if (!is.null(run_before)) {
        if (is.list(run_before)) {
          # select only the functions for missing variables
          run_before <- unique(unlist(run_before[missing_vars]))
        }
        run_before_msg <- paste(
          "You may need to apply some or all of the following functions first:",
          paste(paste0(run_before), collapse = ", ")
        )
      } else {
        run_before_msg <- NULL
      }
      error_msg <- paste(c(error_msg, run_before_msg, advice), collapse = "\n")
      cli::cli_abort(error_msg, call = calling_env)
    } else {
      return(all(var_check))
    }
  } else {
    return(all(var_check))
  }
}

#' Ensure that modeldata contains all necessary variables
#'
#' Throws an error if necessary variables are missing, and replaces optional
#' variables with defaults indicating their absence in the stan model.
#'
#' The reason why defaults are inserted for missing optional variables is that
#' stan does not explicitly allow default variables, so they are specified with
#' a little hack by using zero-length dimensions.
#'
#' @param modeldata A `list` with all data variables (including priors) to be
#'   passed on to stan, alongside additional meta information and descriptions.
#' @param defaults A `list` with default values to be used for modeldata
#'   variables if not supplied in modeldata. For example, `numeric(0)` will
#'   often be supplied for optional parameters.
#' @param data A list with observation data as provided by [sewer_data()].
#' @param assumptions A list with assumptions as provided by
#'   [sewer_assumptions()].
#'
#' @return If no error is thrown due to missing mandatory variables, the same
#'   modeldata object is returned again, where optional variables have been
#'   replaced by zero-length dimension defaults for stan.
modeldata_validate <- function(modeldata,
                               data = list(),
                               assumptions = list(),
                               defaults = modeldata_defaults()) {
  if (!class(modeldata) == "modeldata") {
    cli::cli_abort("Please supply a modeldata object.")
  }

  # check for conflicting data
  for (comp in names(modeldata$.sewer_data)) {
    overlaps_data <- intersect(
      names(data), names(modeldata$.sewer_data[[comp]])
    )
    for (o in overlaps_data) {
      if (!is.null(data[[o]]) &&
        !is.null(modeldata$.sewer_data[[comp]][o]) &&
        !identical(data[[o]], modeldata$.sewer_data[[comp]][o][[1]])) {
        rlang::cli_abort(paste0(
          "You provided different data for `", o,
          "` in sewer_data() and ", comp, "()."
        ))
      }
    }
  }

  # check for conflicting assumptions
  for (comp in names(modeldata$.sewer_assumptions)) {
    overlaps_assumptions <- intersect(
      names(assumptions), names(modeldata$.sewer_assumptions[[comp]])
    )
    for (o in overlaps_assumptions) {
      if (!is.null(assumptions[[o]]) &&
        !is.null(modeldata$.sewer_assumptions[[comp]][o]) &&
        !identical(
          assumptions[[o]], modeldata$.sewer_assumptions[[comp]][o][[1]]
        )) {
        cli::cli_abort(paste0(
          "You provided different assumptions for `", o,
          "` in sewer_assumptions() and ", comp, "()."
        ))
      }
    }
  }

  for (i in 1:2) { # update two times to catch all lazy evals
    modeldata <- modeldata_update(
      modeldata,
      data = data, assumptions = assumptions, throw_error = F
    )
  }

  for (modeldata_sub in list(modeldata$.metainfo, modeldata, modeldata$.init)) {
    for (name in names(modeldata_sub)) {
      if ("tbp" %in% class(modeldata_sub[[name]])) {
        message_data <- c()
        message_assumptions <- c()
        if (length(modeldata_sub[[name]]$required_data) > 0) {
          message_data <- paste(
            "Data (also via sewer_data()):",
            paste(
              "`", modeldata_sub[[name]]$required_data, "`",
              sep = "", collapse = ", "
            )
          )
        }
        if (length(modeldata_sub[[name]]$required_assumptions) > 0) {
          message_assumptions <- paste(
            "Assumptions (also via sewer_assumptions()):",
            paste(
              "`", modeldata_sub[[name]]$required_assumptions, "`",
              sep = "", collapse = ", "
            )
          )
        }
        cli::cli_abort(c(
          paste0(
            "Please provide the following information to ", name, "():"
          ),
          message_data,
          message_assumptions
        ))
      }
    }
  }

  for (modeldata_sub in list(modeldata$.metainfo, modeldata, modeldata$.init)) {
    for (name in names(modeldata_sub)) {
      if (any(c("tbe", "tbc") %in% class(modeldata_sub[[name]]))) {
        error_message <- paste0(
          "There is still information missing to compute ",
          "`", name, "`"
        )
        if (!is.null(modeldata_sub[[name]]$advice)) {
          error_message <- c(error_message, modeldata_sub[[name]]$advice)
        }
        cli::cli_abort(error_message)
      }
    }
  }

  modeldata <- modeldata_update(
    modeldata,
    data = data, assumptions = assumptions
  )

  for (i in seq_along(defaults)) {
    levels <- stringr::str_split(names(defaults)[[i]], "\\$")[[1]]
    modeldata <- default_list_nested(
      modeldata,
      levels = levels, default = defaults[[i]]
    )
  }

  for (check in modeldata$.checks) {
    check(modeldata)
  }

  return(modeldata)
}

modeldata_combine <- function(...) {
  modeldata_sets <- list(...)

  lapply(modeldata_sets, function(md) {
    if (!class(md) == "modeldata") {
      cli::cli_abort("Please supply a modeldata object.")
    }
  })

  modeldata_combined <- do.call(
    c, lapply(modeldata_sets, function(x) {
      list_except(
        x, c(
          ".init", ".metainfo", ".checks", ".str",
          ".sewer_data", ".sewer_assumptions"
        )
      )
    })
  )
  class(modeldata_combined) <- "modeldata"
  modeldata_combined$.init <- do.call(
    c, lapply(modeldata_sets, function(x) x$.init)
  )
  modeldata_combined$.metainfo <- do.call(
    c, lapply(modeldata_sets, function(x) x$.metainfo)
  )
  modeldata_combined$.checks <- do.call(
    c, lapply(modeldata_sets, function(x) x$.checks)
  )
  modeldata_combined$.str <- list(
    measurements = do.call(
      c, lapply(modeldata_sets, function(x) x$.str$measurements),
    ),
    sampling = do.call(
      c, lapply(modeldata_sets, function(x) x$.str$sampling),
    ),
    sewage = do.call(
      c, lapply(modeldata_sets, function(x) x$.str$sewage),
    ),
    shedding = do.call(
      c, lapply(modeldata_sets, function(x) x$.str$shedding),
    ),
    infections = do.call(
      c, lapply(modeldata_sets, function(x) x$.str$infections),
    )
  )
  class(modeldata_combined$.str) <- c("modelstructure")
  modeldata_combined$.sewer_data <- do.call(
    c, lapply(modeldata_sets, function(x) x$.sewer_data)
  )
  modeldata_combined$.sewer_assumptions <- do.call(
    c, lapply(modeldata_sets, function(x) x$.sewer_assumptions)
  )
  modeldata_combined <- modeldata_update(
    modeldata_combined,
    throw_error = FALSE
  )
  return(modeldata_combined)
}

modeldata_update <- function(modeldata,
                             data = list(), assumptions = list(),
                             throw_error = TRUE) {
  if (!class(modeldata) == "modeldata") {
    cli::cli_abort("Please supply a modeldata object.")
  }

  modeldata <- modeldata_update_metainfo(modeldata)

  # meta info
  all_vars <- names(modeldata$.metainfo)
  for (name in all_vars) {
    if ("tbe" %in% class(modeldata$.metainfo[[name]])) {
      modeldata$.metainfo[[name]] <- solve(
        modeldata$.metainfo[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbc" %in% class(modeldata$.metainfo[[name]])) {
      modeldata <- solve(
        modeldata$.metainfo[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbp" %in% class(modeldata$.metainfo[[name]])) {
      modeldata <- solve(
        modeldata$.metainfo[[name]],
        modeldata = modeldata,
        data = data,
        assumptions = assumptions
      )
    }
  }

  # main modeldata
  all_vars <- setdiff(names(modeldata), c(".init", ".metainfo", ".checks"))
  for (name in all_vars) {
    if ("tbe" %in% class(modeldata[[name]])) {
      modeldata[[name]] <- solve(
        modeldata[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbc" %in% class(modeldata[[name]])) {
      modeldata <- solve(
        modeldata[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbp" %in% class(modeldata[[name]])) {
      modeldata <- solve(
        modeldata[[name]],
        modeldata = modeldata,
        data = data,
        assumptions = assumptions
      )
    }
  }

  # init
  all_vars <- names(modeldata$.init)
  for (name in all_vars) {
    if ("tbe" %in% class(modeldata$.init[[name]])) {
      modeldata$.init[[name]] <- solve(
        modeldata$.init[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbc" %in% class(modeldata$.init[[name]])) {
      modeldata <- solve(
        modeldata$.init[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbp" %in% class(modeldata$.init[[name]])) {
      modeldata <- solve(
        modeldata$.init[[name]],
        modeldata = modeldata,
        data = data,
        assumptions = assumptions
      )
    }
  }

  modeldata <- modeldata_update_metainfo(modeldata)

  return(modeldata)
}



#' Print a `modeldata` object
#'
#' @param type If "structure", the model structure is printed. If "data", the
#'   modeldata content (data, inits, metainfo) is printed.
#' @export
print.modeldata <- function(x, type = "structure") {
  if (type == "structure") {
    print(x$.str)
  } else if (type == "data") {
    print.default(x[!names(x) %in% c(
      ".str", ".checks", ".sewer_data", ".sewer_assumptions"
      )])
  } else {
    print.default(x)
  }
}

#' Print a `modelstructure` object
#'
#' @export
print.modelstructure <- function(.str) {
  output <- sapply(names(.str), function(module) {
    if (length(.str[[module]]) > 0) {
      comp_output <- sapply(names(.str[[module]]), function(component) {
        if (length(.str[[module]][[component]][[1]]) > 0) {
          a <- .str[[module]][[component]][[1]]
          details <- paste0(
            " (", paste(paste0(names(a), " = ", a), collapse = ", "), ")"
            )
          paste0(" |- ", names(.str[[module]][[component]])[1], details)
        } else {
          paste0(" |- ", names(.str[[module]][[component]])[1])
        }
      })
      paste0(module, "\n", paste(comp_output, collapse = "\n"))
    } else {
      return(NULL)
    }
  })
  cat(paste(output[!sapply(output,is.null)], collapse = "\n\n"))
}
