## ---------------------------------------------------------------
##                           Stan utils                          -
## ---------------------------------------------------------------

#' Title
#'
#' @param model_filename
#' @param standata
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
    standata = NULL,
    model_folder = "stan",
    profile = TRUE,
    threads = FALSE,
    force_recompile = FALSE) {
  model_def <- list()

  if (is.null(model_filename)) {
    if (is.null(standata)) {
      abort(
        c(
          "Please either provide ",
          "`model_filename` and `model_folder` (i.e. path to a stan model)",
          "`standata` (so that a suitable stan model can be inferred)"
        )
      )
    }
    if (standata$meta_info$R_estimate_approach == "splines") {
      model_filename <- "wastewater_Re_splines.stan"
    } else if (standata$meta_info$R_estimate_approach == "ets") {
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
set_fit_opts <- function(sampler_opts = sampler_opts(), fitted = TRUE) {
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
set_sampler_opts <- function(
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
    show_messages = T,
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

#' Check standata object for required variables
#'
#' `standata_check` accepts various additional arguments which are turned into a
#' sophisticated error message about missing variables
#'
#' @inheritParams standata_validate
#' @param required Character vector of required variables
#' @param descriptions Optional, a list with names corresponding to all or some
#'   of the required variables, and values corresponding to additional
#'   descriptions of the variables
#' @param run_before Optional, either a character vector, or a `list` with names
#'   corresponding to all of the required variables, and values
#'   corresponding to functions that will add these variables to standata.
#' @param advice Additional message with advice on how to solve error.
#' @param throw_error Should an error be thrown (with verbatim error message),
#'   or should `TRUE`/`FALSE` be returned. Default is `TRUE`.
#'
#' @return If throw_error is `TRUE`, nothing is returned, otherwise a logical is
#'   returned (`TRUE` if all required variables are present)
standata_check <- function(standata,
                           required,
                           descriptions = standata_descriptions(),
                           run_before = standata_var_requirements(),
                           advice = NULL,
                           throw_error = T,
                           calling_env = rlang::caller_env()) {
  var_check <- check_list_nested(standata, required)
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
        "but currently not present in standata:\n\n",
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

#' Ensure that standata contains all necessary variables
#'
#' Throws an error if necessary variables are missing, and replaces optional
#' variables with defaults indicating their absence in the stan model.
#'
#' The reason why defaults are inserted for missing optional variables is that
#' stan does not explicitly allow default variables, so they are specified with
#' a little hack by using zero-length dimensions.
#'
#' @param standata A `list` with all data variables (including priors) to be
#'   passed on to stan, alongside additional meta information and descriptions.
#' @param model_def A `list` with information about the stan model to be fitted
#'   and a function that returns the CmdStanModel object.
#' @param defaults A `list` with default values to be used for standata
#'   variables if not supplied in standata. For example, `numeric(0)` will often
#'   be supplied for optional parameters.
#'
#' @return If no error is thrown due to missing mandatory variables, the same
#'   standata object is returned again, where optional variables have been
#'   replaced by zero-length dimension defaults for stan.
standata_validate <- function(standata, model_def,
                              defaults = standata_defaults()) {
  standata <- standata_update(standata)

  for (standata_sub in list(standata$meta_info, standata, standata$init)) {
    for (name in names(standata_sub)) {
      if (any(c("tbe", "tbef") %in% class(standata_sub[[name]]))) {
        abort(paste0(
          "There is still information missing to compute ",
          "`", name, "`"
        ))
      }
    }
  }

  for (i in seq_along(defaults)) {
    levels <- stringr::str_split(names(defaults)[[i]], "\\$")[[1]]
    standata <- default_list_nested(
      standata,
      levels = levels, default = defaults[[i]]
    )
  }

  return(standata)
}

standata_combine <- function(...) {
  standata_sets <- list(...)
  standata_combined <- do.call(
    c, lapply(standata_sets, function(x) {
      list_except(x, c("init", "meta_info"))
    })
  )
  standata_combined$init <- do.call(
    c, lapply(standata_sets, function(x) x$init)
  )
  standata_combined$meta_info <- do.call(
    c, lapply(standata_sets, function(x) x$meta_info)
  )
  standata_combined <- standata_update(standata_combined, throw_error = FALSE)
  return(standata_combined)
}

standata_update <- function(standata, throw_error = TRUE) {
  standata <- standata_update_metainfo(standata)

  # meta info
  all_vars <- names(standata$meta_info)
  for (name in all_vars) {
    if ("tbe" %in% class(standata$meta_info[[name]])) {
      standata$meta_info[[name]] <- solve(
        standata$meta_info[[name]],
        standata = standata, throw_error = throw_error
      )
    } else if ("tbef" %in% class(standata$meta_info[[name]])) {
      standata <- solve(
        standata$meta_info[[name]],
        standata = standata, throw_error = throw_error
      )
    }
  }

  # main standata
  all_vars <- setdiff(names(standata), c("init", "meta_info"))
  for (name in all_vars) {
    if ("tbe" %in% class(standata[[name]])) {
      standata[[name]] <- solve(
        standata[[name]],
        standata = standata, throw_error = throw_error
      )
    } else if ("tbef" %in% class(standata[[name]])) {
      standata <- solve(
        standata[[name]],
        standata = standata, throw_error = throw_error
      )
    }
  }

  # init
  all_vars <- names(standata$init)
  for (name in all_vars) {
    if ("tbe" %in% class(standata$init[[name]])) {
      standata$init[[name]] <- solve(
        standata$init[[name]],
        standata = standata, throw_error = throw_error
      )
    } else if ("tbef" %in% class(standata$init[[name]])) {
      standata <- solve(
        standata$init[[name]],
        standata = standata, throw_error = throw_error
      )
    }
  }

  return(standata)
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
stan_prior <- function(param, dist = "normal", ...) {
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

tbe <- function(r_expr, required = c(), calling_env = rlang::caller_env()) {
  if (standata_check(calling_env$standata, required, throw_error = F)) {
    return(r_expr)
  } else {
    lazy_r <- lazyeval::lazy(r_expr)
    # rm(ets_diff, envir = calling_env)
    lazy_r$env <- calling_env
    tbe_o <- list()
    tbe_o$lazy_r <- lazy_r
    tbe_o$calling_f <- tail(rlang::trace_back()$call, 2)[[1]]
    tbe_o$required <- required
    class(tbe_o) <- "tbe"
    return(tbe_o)
  }
}

print.tbe <- function(x) {
  print(paste("Waiting for meta-information:", expr_text(x$lazy_r$expr)))
}

solve.tbe <- function(x, standata, throw_error = TRUE) {
  if (standata_check(
    standata, x$required,
    calling_env = x$calling_f, throw_error = throw_error
  )) {
    return(lazyeval::lazy_eval(x$lazy_r, list(standata = standata)))
  } else {
    return(x)
  }
}

tbef <- function(f_name, f_expr, required = c(),
                 calling_env = rlang::caller_env()) {
  calling_standata <- calling_env$standata

  f_lazy <- lazyeval::lazy(f_expr)
  f_func <- function(standata) {
    f_lazy$env$standata <- standata
    lazyeval::lazy_eval(f_lazy)
    return(f_lazy$env$standata)
  }
  rm(standata, envir = calling_env)
  calling_env$f_lazy <- f_lazy
  environment(f_func) <- calling_env
  if (standata_check(calling_standata, required = required, throw_error = F)) {
    new_standata <- f_func(calling_standata)
  } else {
    tbef_o <- list(name = f_name, func = f_func)
    tbef_o$calling_f <- tail(rlang::trace_back()$call, 2)[[1]]
    tbef_o$required <- required
    class(tbef_o) <- "tbef"
    new_standata <- calling_standata
    new_standata[[f_name]] <- tbef_o
  }
  return(new_standata)
}

print.tbef <- function(x) {
  print(paste(
    "Waiting for the following meta-information to compute values:",
    paste(x$required, collapse = ", ")
  ))
}

solve.tbef <- function(x, standata, throw_error = TRUE) {
  if (standata_check(
    standata, x$required,
    calling_env = x$calling_f,
    throw_error = throw_error
  )) {
    standata <- x$func(standata)
    standata[[x$name]] <- NULL
  }
  return(standata)
}
