#' Template function
#'
#' @description This function does not actually do anything. It only serves as a
#'   template for the documentation of other functions using inheritance of
#'   params.
#'
#' @param modeldata A `modeldata` object to which the above model specifications
#'   should be added. Default is an empty model given by [modeldata_init()]. Can
#'   also be an already partly specified model returned by other `EpiSewer`
#'   model helpers.
#'
#' @return Nothing
template_model_helpers <- function(modeldata) { }

#' Construct an unspecified EpiSewer model
#'
#' @return A `modeldata` object containing data and specifications of the model
#'   to be fitted. Can be passed on to other `EpiSewer` model helpers to add
#'   further data and model specifications.
#'
#' @return The `modeldata` object also includes information about parameter
#'   initialization (`init`), meta data (`meta_info`), and checks to be
#'   performed before model fitting (`checks`).
#' @export
#'
#' @examples
#' modeldata_init()
modeldata_init <- function() {
  modeldata <- list()
  modeldata$init <- list()
  modeldata$meta_info <- list()
  modeldata$checks <- list()
  class(modeldata) <- "modeldata"
  return(modeldata)
}

#' Verification of modeldata objects
verify_is_modeldata <- function(modeldata, arg_name) {
  if (!all(c("init", "meta_info", "checks") %in% names(modeldata))) {
    error_msg <- paste0(
      "The argument `", arg_name, "` expects EpiSewer model data as input. ",
      "Please use a function of the form `", arg_name, "_`", " to specify ",
      "this argument."
    )
    all_fs <- names(rlang::ns_env("EpiSewer"))
    arg_fs <- all_fs[stringr::str_detect(all_fs, paste0("^", arg_name, "_"))]
    if (length(arg_fs)>0) {
      error_msg <- paste(error_msg, "Available functions:")
      error_msg <- c(error_msg, paste0(arg_fs, "()"))
    }
    rlang::abort(error_msg, call = rlang::caller_env())
  }
}

#' Get model helpers for a component (internal function)
#'
#' Get all available model helpers for a given component. Argument defaults are
#' such that helpers are returned as a comma-separated list of roxygen function
#' links.
#'
#' @param component A character vector with the name of the component
#' @param collapse A character used to collapse helpers into one string.
#' Default is NULL (do not collapse).
#' @param prefix A prefix to add to each function name
#' @param suffix A suffix to add to each function name
#'
#' @return A character vector with all available model helper functions for
#' the component. If collapse is not NULL, the helper functions are collapsed
#' into a single string.
component_helpers_ <- function(component, collapse = "\n",
                               prefix = "- [", suffix = "()]") {
  if (!component %in% all_components()) {
    rlang::abort(c(paste("No valid component provided.",
                         "Must be one out of:"),
                   all_components()))
  }
  all_fs <- names(rlang::ns_env("EpiSewer"))
  arg_fs <- all_fs[stringr::str_detect(all_fs, paste0("^", component, "_"))]
  if (length(arg_fs)>0) {
    helpers <- paste0(prefix, arg_fs, suffix)
  } else {
    helpers <- c()
  }
  if (!is.null(collapse)) {
    return(paste(helpers, collapse = collapse))
  } else {
    return(helpers)
  }
}

#' Get model helpers for a component
#'
#' @inheritParams component_helpers_
#'
#' @return A character vector with all available model helper functions for
#' the component.
#' @export
#'
#' @examples
#' component_helpers("R")
component_helpers <- function(component) {
  return(component_helpers_(
    component, collapse = NULL, prefix = "", suffix = "()"
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

#' Check modeldata object for required variables
#'
#' `modeldata_check` accepts various additional arguments which are turned into a
#' sophisticated error message about missing variables
#'
#' @inheritParams modeldata_validate
#' @param required Character vector of required variables
#' @param descriptions Optional, a list with names corresponding to all or some
#'   of the required variables, and values corresponding to additional
#'   descriptions of the variables
#' @param run_before Optional, either a character vector, or a `list` with names
#'   corresponding to all of the required variables, and values
#'   corresponding to functions that will add these variables to modeldata.
#' @param advice Additional message with advice on how to solve error.
#' @param throw_error Should an error be thrown (with verbatim error message),
#'   or should `TRUE`/`FALSE` be returned. Default is `TRUE`.
#'
#' @return If throw_error is `TRUE`, nothing is returned, otherwise a logical is
#'   returned (`TRUE` if all required variables are present)
modeldata_check <- function(modeldata,
                            required,
                            descriptions = modeldata_descriptions(),
                            run_before = modeldata_var_requirements(),
                            advice = NULL,
                            throw_error = TRUE,
                            calling_env = rlang::caller_env()) {

  if (!class(modeldata) == "modeldata") {
    rlang::abort("Please supply a modeldata object.")
  }

  var_check <- check_list_nested(modeldata, required)
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
      rlang::abort(error_msg, call = calling_env)
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
#' @param model_def A `list` with information about the stan model to be fitted
#'   and a function that returns the CmdStanModel object.
#' @param defaults A `list` with default values to be used for modeldata
#'   variables if not supplied in modeldata. For example, `numeric(0)` will often
#'   be supplied for optional parameters.
#' @param data
#' @param assumptions
#'
#' @return If no error is thrown due to missing mandatory variables, the same
#'   modeldata object is returned again, where optional variables have been
#'   replaced by zero-length dimension defaults for stan.
modeldata_validate <- function(modeldata,
                               data = list(),
                               assumptions = list(),
                               model_def,
                               defaults = modeldata_defaults()) {
  if (!class(modeldata) == "modeldata") {
    rlang::abort("Please supply a modeldata object.")
  }

  # check for conflicting data
  for (comp in names(modeldata$sewer_data)) {
    overlaps_data <- intersect(
      names(data), names(modeldata$sewer_data[[comp]])
    )
    for (o in overlaps_data) {
      if (!is.null(data[[o]]) &&
          !is.null(modeldata$sewer_data[[comp]][o]) &&
          !identical(data[[o]], modeldata$sewer_data[[comp]][o][[1]])) {
        rlang::cli_abort(paste0(
          "You provided different data for `", o,
          "` in sewer_data() and ", comp, "()."
        ))
      }
    }
  }

  # check for conflicting assumptions
  for (comp in names(modeldata$sewer_assumptions)) {
    overlaps_assumptions <- intersect(
      names(assumptions), names(modeldata$sewer_assumptions[[comp]])
    )
    for (o in overlaps_assumptions) {
      if (!is.null(assumptions[[o]]) &&
          !is.null(modeldata$sewer_assumptions[[comp]][o]) &&
          !identical(assumptions[[o]], modeldata$sewer_assumptions[[comp]][o][[1]])) {
        rlang::abort(paste0(
          "You provided different assumptions for `", o,
          "` in sewer_assumptions() and ", comp, "()."
        ))
      }
    }
  }

  for (i in 1:2) { # update two times to catch all lazy evals
    modeldata <- modeldata_update(
      modeldata, data = data, assumptions = assumptions, throw_error = F
    )
  }

  for (modeldata_sub in list(modeldata$meta_info, modeldata, modeldata$init)) {
    for (name in names(modeldata_sub)) {
      if ("tbp" %in% class(modeldata_sub[[name]])) {
        message_data <- c()
        message_assumptions <- c()
        if (length(modeldata_sub[[name]]$required_data) > 0) {
          message_data <- paste(
            "Data (also via sewer_data()):",
            paste("`", modeldata_sub[[name]]$required_data, "`", sep = "", collapse = ", ")
            )
        }
        if (length(modeldata_sub[[name]]$required_assumptions) > 0) {
          message_assumptions <- paste(
            "Assumptions (also via sewer_assumptions()):",
            paste("`", modeldata_sub[[name]]$required_assumptions, "`", sep = "", collapse = ", ")
            )
        }
        rlang::abort(c(paste0(
          "Please provide the following information to ", name, "():"
        ),
        message_data,
        message_assumptions
        ))
      }
    }
  }

  for (modeldata_sub in list(modeldata$meta_info, modeldata, modeldata$init)) {
    for (name in names(modeldata_sub)) {
      if (any(c("tbe", "tbc") %in% class(modeldata_sub[[name]]))) {
        rlang::abort(paste0(
          "There is still information missing to compute ",
          "`", name, "`"
        ))
      }
    }
  }

  modeldata <- modeldata_update(
    modeldata, data = data, assumptions = assumptions
    )

  for (i in seq_along(defaults)) {
    levels <- stringr::str_split(names(defaults)[[i]], "\\$")[[1]]
    modeldata <- default_list_nested(
      modeldata,
      levels = levels, default = defaults[[i]]
    )
  }

  for (check in modeldata$checks) {
    check(modeldata)
  }

  return(modeldata)
}

modeldata_combine <- function(...) {
  modeldata_sets <- list(...)

  lapply(modeldata_sets, function(md) {
    if (!class(md) == "modeldata") {
      rlang::abort("Please supply a modeldata object.")
    }
  })

  modeldata_combined <- do.call(
    c, lapply(modeldata_sets, function(x) {
      list_except(
        x, c("init", "meta_info", "checks", "sewer_data", "sewer_assumptions")
        )
    })
  )
  class(modeldata_combined) <- "modeldata"
  modeldata_combined$init <- do.call(
    c, lapply(modeldata_sets, function(x) x$init)
  )
  modeldata_combined$meta_info <- do.call(
    c, lapply(modeldata_sets, function(x) x$meta_info)
  )
  modeldata_combined$checks <- do.call(
    c, lapply(modeldata_sets, function(x) x$checks)
  )
  modeldata_combined$sewer_data <- do.call(
    c, lapply(modeldata_sets, function(x) x$sewer_data)
  )
  modeldata_combined$sewer_assumptions <- do.call(
    c, lapply(modeldata_sets, function(x) x$sewer_assumptions)
  )
  modeldata_combined <- modeldata_update(
    modeldata_combined, throw_error = FALSE
    )
  return(modeldata_combined)
}

modeldata_update <- function(modeldata,
                             data = list(), assumptions = list(),
                             throw_error = TRUE) {
  if (!class(modeldata) == "modeldata") {
    rlang::abort("Please supply a modeldata object.")
  }

  modeldata <- modeldata_update_metainfo(modeldata)

  # meta info
  all_vars <- names(modeldata$meta_info)
  for (name in all_vars) {
    if ("tbe" %in% class(modeldata$meta_info[[name]])) {
      modeldata$meta_info[[name]] <- solve(
        modeldata$meta_info[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbc" %in% class(modeldata$meta_info[[name]])) {
      modeldata <- solve(
        modeldata$meta_info[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbp" %in% class(modeldata$meta_info[[name]])) {
      modeldata <- solve(
        modeldata$meta_info[[name]],
        modeldata = modeldata,
        data = data,
        assumptions = assumptions
      )
    }
  }

  # main modeldata
  all_vars <- setdiff(names(modeldata), c("init", "meta_info", "checks"))
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
  all_vars <- names(modeldata$init)
  for (name in all_vars) {
    if ("tbe" %in% class(modeldata$init[[name]])) {
      modeldata$init[[name]] <- solve(
        modeldata$init[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbc" %in% class(modeldata$init[[name]])) {
      modeldata <- solve(
        modeldata$init[[name]],
        modeldata = modeldata, throw_error = throw_error
      )
    } else if ("tbp" %in% class(modeldata$init[[name]])) {
      modeldata <- solve(
        modeldata$init[[name]],
        modeldata = modeldata,
        data = data,
        assumptions = assumptions
      )
    }
  }

  modeldata <- modeldata_update_metainfo(modeldata)

  return(modeldata)
}
