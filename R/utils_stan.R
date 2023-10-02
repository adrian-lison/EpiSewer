## ---------------------------------------------------------------
##                           Stan utils                          -
## ---------------------------------------------------------------

#' Title
#'
#' @param model_filename
#' @param modeldata
#' @param model_folder
#' @param profile
#' @param threads
#' @param force_recompile
#'
#' @return
#' @export
#'
#' @examples
get_stan_model <- function(
    model_filename = NULL,
    modeldata = NULL,
    model_folder = "stan",
    profile = TRUE,
    threads = FALSE,
    force_recompile = FALSE) {
  model_def <- list()

  if (is.null(model_filename)) {
    if (is.null(modeldata)) {
      abort(
        c(
          "Please either provide ",
          "`model_filename` and `model_folder` (i.e. path to a stan model)",
          "or `modeldata` (so that a suitable stan model can be inferred)."
        )
      )
    }
    if (modeldata$meta_info$R_estimate_approach == "splines") {
      model_filename <- "wastewater_Re_splines.stan"
    } else if (modeldata$meta_info$R_estimate_approach == "ets") {
      model_filename <- "wastewater_Re.stan"
    } else if (modeldata$meta_info$R_estimate_approach == "rw") {
      model_filename <- "wastewater_Re.stan"
    } else {
      abort(
        paste(
          "No suitable model available,",
          "please supply a stan model yourself",
          "using `get_stan_model`, or open an issue."
        )
      )
    }
  }

  model_def[["model_filename"]] <- model_filename
  model_def[["model_folder"]] <- model_folder
  model_def[["force_recompile"]] <- force_recompile
  model_def[["threads"]] <- threads
  model_def[["profile"]] <- profile

  model_def <- update_compiled_stanmodel(model_def, force_recompile)

  return(model_def)
}

#' Title
#'
#' @param sampler_opts
#' @param fitted
#'
#' @return
#' @export
#'
#' @examples
set_fit_opts <- function(sampler = sampler_stan_mcmc(), fitted = TRUE) {
  opts <- as.list(environment())
  return(opts)
}

#' Title
#'
#' @param chains
#' @param iter_warmup
#' @param iter_sampling
#' @param adapt_delta
#' @param max_treedepth
#' @param step_size
#' @param parallel_chains
#' @param threads_per_chain
#' @param seed
#' @param refresh
#' @param show_messages
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sampler_stan_mcmc <- function(
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    max_treedepth = 15,
    step_size = 0.1,
    parallel_chains = NULL,
    threads_per_chain = 1,
    seed = 0,
    refresh = 200,
    show_messages = TRUE,
    show_exceptions = FALSE,
    ...) {
  opts <- c(as.list(environment()), list(...))
  if (opts$threads_per_chain == 1) {
    opts$threads_per_chain <- NULL
  }
  return(opts)
}

update_compiled_stanmodel <- function(model_def, force_recompile = FALSE) {
  n_models <- length(model_def$model_filename)

  model_path_list <- lapply(1:n_models, function(i) {
    system.file(model_def$model_folder,
                model_def$model_filename[[i]], package = "EpiSewer")
  })
  include_paths_list <- lapply(1:n_models, function(i) {
    system.file(model_def$model_folder, package = "EpiSewer")
  }) # identical

  if (!model_def$profile) {
    stan_no_profiling_list <- lapply(1:n_models, function(i) {
      write_stan_files_no_profiling(
        model_path_list[[i]],
        include_paths_list[[i]]
      )
    })
    model_path_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$model
    })
    include_paths_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$include_paths
    })
  }

  cpp_options <- list()
  if (model_def$threads) {
    cpp_options[["stan_threads"]] <- TRUE
  }

  stanmodel_list <- lapply(1:n_models, function(i) {
    try(cmdstanr::cmdstan_model(model_path_list[[i]],
      include_paths = include_paths_list[[i]],
      dir = dirname(model_path_list[[i]]),
      cpp_options = cpp_options,
      force_recompile = force_recompile
    ))
  })

  model_def[["get_stan_model"]] <- lapply(1:n_models, function(i) {
    function() {
      return(stanmodel_list[[i]])
    }
  })

  return(model_def)
}

verify_is_modeldata <- function(modeldata, arg_name) {
  if (!all(c("init", "meta_info", "checks") %in% names(modeldata))) {
    error_msg <- paste0(
      "The argument `", arg_name, "` expects EpiSewer model data as input. ",
      "Please use a function of the form `", arg_name, "_`", " to specify ",
      "this argument."
    )
    all_fs <- names(ns_env("EpiSewer"))
    arg_fs <- all_fs[stringr::str_detect(all_fs, paste0("^", arg_name, "_"))]
    if (length(arg_fs)>0) {
      error_msg <- paste(error_msg, "Available functions:")
      error_msg <- c(error_msg, paste0(arg_fs, "()"))
    }
    abort(error_msg, call = caller_env())
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
    abort(c(paste("No valid component provided.",
                "Must be one out of:"),
                all_components()))
  }
  all_fs <- names(ns_env("EpiSewer"))
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
    abort("Please supply a modeldata object.")
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
      abort(error_msg, call = calling_env)
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
#'
#' @return If no error is thrown due to missing mandatory variables, the same
#'   modeldata object is returned again, where optional variables have been
#'   replaced by zero-length dimension defaults for stan.
modeldata_validate <- function(modeldata, model_def,
                              defaults = modeldata_defaults()) {
  if (!class(modeldata) == "modeldata") {
    abort("Please supply a modeldata object.")
  }

  modeldata <- modeldata_update(modeldata)

  for (modeldata_sub in list(modeldata$meta_info, modeldata, modeldata$init)) {
    for (name in names(modeldata_sub)) {
      if (any(c("tbe", "tbc") %in% class(modeldata_sub[[name]]))) {
        abort(paste0(
          "There is still information missing to compute ",
          "`", name, "`"
        ))
      }
    }
  }

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
      abort("Please supply a modeldata object.")
    }
  })

  modeldata_combined <- do.call(
    c, lapply(modeldata_sets, function(x) {
      list_except(x, c("init", "meta_info", "checks"))
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
  modeldata_combined <- modeldata_update(modeldata_combined, throw_error = FALSE)
  return(modeldata_combined)
}

modeldata_update <- function(modeldata, throw_error = TRUE) {
  if (!class(modeldata) == "modeldata") {
    abort("Please supply a modeldata object.")
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
    }
  }

  return(modeldata)
}

#' Title
#'
#' @param param
#' @param dist
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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


#' To be assumed
tba <- function(r_expr, required = c(), calling_env = rlang::caller_env()) {
  if (required %in% calling_env$assumptions) {
    return(r_expr)
  } else {
    lazy_r <- lazyeval::lazy(r_expr)
    lazy_r$env <- calling_env
    tba_o <- list()
    tba_o$lazy_r <- lazy_r
    tba_o$calling_f <- tail(rlang::trace_back()$call, 2)[[1]]
    tba_o$required <- required
    class(tba_o) <- "tba"
    return(tba_o)
  }
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
print.tba <- function(x) {
  print(paste("Requires assumption:", expr_text(x$lazy_r$expr)))
}

solve.tba <- function(x, assumptions) {
  if (x$required %in% assumptions) {
    return(lazyeval::lazy_eval(x$lazy_r, list(assumptions = assumptions)))
  } else {
    return(x)
  }
}

#' To be evaluated
tbe <- function(r_expr, required = c(), calling_env = rlang::caller_env()) {
  if (modeldata_check(calling_env$modeldata, required, throw_error = FALSE)) {
    return(r_expr)
  } else {
    lazy_r <- lazyeval::lazy(r_expr)
    lazy_r$env <- calling_env
    tbe_o <- list()
    tbe_o$lazy_r <- lazy_r
    tbe_o$calling_f <- tail(rlang::trace_back()$call, 2)[[1]]
    tbe_o$required <- required
    class(tbe_o) <- "tbe"
    return(tbe_o)
  }
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
print.tbe <- function(x) {
  print(paste("Waiting for information:", expr_text(x$lazy_r$expr)))
}

solve.tbe <- function(x, modeldata, throw_error = TRUE) {
  if (modeldata_check(
    modeldata, x$required,
    calling_env = x$calling_f, throw_error = throw_error
  )) {
    return(lazyeval::lazy_eval(x$lazy_r, list(modeldata = modeldata)))
  } else {
    return(x)
  }
}

#' To be computed
tbc <- function(f_name, f_expr, required = c(),
                 calling_env = rlang::caller_env()) {
  calling_modeldata <- calling_env$modeldata

  f_lazy <- lazyeval::lazy(f_expr)
  f_func <- function(modeldata) {
    f_lazy$env$modeldata <- modeldata
    lazyeval::lazy_eval(f_lazy)
    return(f_lazy$env$modeldata)
  }
  rm(modeldata, envir = calling_env)
  calling_env$f_lazy <- f_lazy
  environment(f_func) <- calling_env
  if (modeldata_check(
    calling_modeldata, required = required, throw_error = FALSE)
    ) {
    new_modeldata <- f_func(calling_modeldata)
  } else {
    tbc_o <- list(name = f_name, func = f_func)
    tbc_o$calling_f <- tail(rlang::trace_back()$call, 2)[[1]]
    tbc_o$required <- required
    class(tbc_o) <- "tbc"
    new_modeldata <- calling_modeldata
    new_modeldata[[f_name]] <- tbc_o
  }
  return(new_modeldata)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
print.tbc <- function(x) {
  print(paste(
    "Waiting for the following information to compute values:",
    paste(x$required, collapse = ", ")
  ))
}

solve.tbc <- function(x, modeldata, throw_error = TRUE) {
  if (modeldata_check(
    modeldata, x$required,
    calling_env = x$calling_f,
    throw_error = throw_error
  )) {
    modeldata <- x$func(modeldata)
    modeldata[[x$name]] <- NULL
  }
  return(modeldata)
}
