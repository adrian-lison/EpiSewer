#' Check if all required data and assumptions have been provided.
#'
#' @param f_env The environment to be checked
#' @param data A list with data
#' @param assumptions A list with assumptions
#' @param required_data A character vector with the names of required data.
#'   Different options for one item can be provided using "|" as a separator.
#'   EpiSewer will then check if at least one of the options is present.
#' @param required_assumptions A character vector with the names of required
#'   assumptions. Different options for one item can be provided using "|" as a
#'   separator. EpiSewer will then check if at least one of the options is
#'   present.
#' @return TRUE if all required data and assumptions are present in `f_env`.
#' @keywords internal
check_provided <- function(f_env = list(), data = list(), assumptions = list(),
                           required_data = c(), required_assumptions = c()) {
  f_env <- update_env(f_env, data, assumptions, required_data, required_assumptions)
  present <- sapply(c(required_data, required_assumptions), function(item) {
    options <- stringr::str_split(string = item, "\\|")[[1]]
    for (o in options) {
      if (o %in% names(f_env) && !is.null(f_env[[o]])) {
        return(TRUE)
      }
    }
    return(FALSE)
  })
  return(all(present))
}

#' Add required data and assumptions to an environment
#'
#' @param f_env The environment to be updated
#' @inheritParams check_provided
#'
#' @details Note that data and assumptions are only assigned to the environment
#'   if the corresponding attributes are not yet present. If already present,
#'   they are not overwritten.
#'
#' @return The updated environment, where all required data and assumptions that
#'   were provided in `data` and `assumptions` were inserted if missing.
#' @keywords internal
update_env <- function(f_env = list(), data = list(), assumptions = list(),
                       required_data = c(), required_assumptions = c()) {
  for (req in required_data) {
    for (r in stringr::str_split(string = req, "\\|")[[1]]) {
      if (!(r %in% names(f_env) && !is.null(f_env[[r]]))) {
        f_env[[r]] <- data[[r]]
      }
    }
  }
  for (req in required_assumptions) {
    for (r in stringr::str_split(string = req, "\\|")[[1]]) {
      if (!(r %in% names(f_env) && !is.null(f_env[[r]]))) {
        f_env[[r]] <- assumptions[[r]]
      }
    }
  }
  return(f_env)
}

#' Copy an environment
#'
#' @param env The environment to be copied
#' @param exclude Character vector with names to be excluded in the copy
#'
#' @details The new environment has the calling environment of `copy_env` as its
#'   parent. This ensure that all functions in the namespace are accessible.
#'
#' @return A copy of the environment
#' @keywords internal
copy_env <- function(env, exclude = c()) {
  copy <- new.env()
  for(name in ls(env, all.names=TRUE)) {
    if (!(name %in% exclude)) {
      assign(name, get(name, env), copy)
    }
  }
  return(copy)
}

#' Lazy-execute a function once all data and assumptions are provided
#'
#' @param f_name Name of the "functions" that will be executed once everything
#'   is provided.
#' @param f_expr Expression with arbitrary R code in which attributes are
#'   assigned to `modeldata.`
#' @param required_data A character vector with the names of required data.
#'   Different options for one item can be provided using "|" as a separator.
#'   EpiSewer will then check if at least one of the options is present.
#' @param required_assumptions A character vector with the names of required
#'   assumptions. Different options for one item can be provided using "|" as a
#'   separator. EpiSewer will then check if at least one of the options is
#'   present.
#' @param modeldata Modeldata object in the calling function.
#' @param calling_env Calling environment, should be `rlang::caller_env()`.
#'
#' @details The difference between `tbp` (to be provided) and `tbc` (to be
#'   computed) is that `tbp` is for promises based on data and assumptions, and
#'   `tbc` is for promises based on internal modeldata variables (typically in
#'   `.metainfo`). This means that there can also be cascades: once a `tbp` is
#'   resolved, this may lead to an update of some metainformation, such that a
#'   `tbc` can be resolved in the next update.
#' @keywords internal
tbp <- function(f_name, f_expr,
                required_data = c(), required_assumptions = c(),
                modeldata = NULL,
                calling_env = rlang::caller_env()) {
  if (is.null(modeldata)) {
    modeldata_temp <- calling_env$modeldata
  } else {
    modeldata_temp <- modeldata
  }

  f_lazy <- lazyeval::lazy(f_expr)
  d_lazy <- as.list(calling_env)
  d_lazy$modeldata <- NULL

  f_func <- function(modeldata,
                     data = list(), assumptions = list(),
                     required_data = c(), required_assumptions = c()) {
    for(n in names(d_lazy)) {
      f_lazy$env <- rlang::env_clone(f_lazy$env)
      assign(n, d_lazy[[n]], f_lazy$env)
    }
    f_lazy$env <- update_env(
      f_lazy$env, data = data, assumptions = assumptions,
      required_data = required_data, required_assumptions = required_assumptions
      )
    f_lazy$env$modeldata <- modeldata

    # register used data
    f_sewer_data <- list()
    f_sewer_data[[f_name]] <- as.list(f_lazy$env)[purrr::list_c(stringr::str_split(required_data, "\\|"))]
    f_lazy$env$modeldata$.sewer_data <- c(
      f_lazy$env$modeldata$.sewer_data, f_sewer_data
    )

    # register used assumptions
    f_sewer_assumptions <- list()
    f_sewer_assumptions[[f_name]] <- as.list(f_lazy$env)[purrr::list_c(stringr::str_split(required_assumptions, "\\|"))]
    f_lazy$env$modeldata$.sewer_assumptions <- c(
      f_lazy$env$modeldata$.sewer_assumptions, f_sewer_assumptions
    )

    # apply function to modeldata
    lazyeval::lazy_eval(f_lazy)

    return(f_lazy$env$modeldata)
  }

  environment(f_func) <- list2env(list(
    f_name = f_name,
    f_lazy = f_lazy,
    d_lazy = d_lazy,
    update_env = update_env,
    required_data = required_data,
    required_assumptions = required_assumptions
  ), copy_env(calling_env, exclude = "modeldata"))

  all_provided <- check_provided(environment(f_func),
    data = list(), assumptions = list(),
    required_data = required_data,
    required_assumptions = required_assumptions
  )

  if (all_provided) {
    modeldata_temp <- f_func(
      modeldata_temp, data = list(), assumptions = list(),
      required_data = required_data, required_assumptions = required_assumptions
    )
  } else {
    tbp_object <- list(
      f_name = f_name,
      f_func = f_func,
      check_provided = check_provided,
      required_data = required_data,
      required_assumptions = required_assumptions,
      calling_f = tail(rlang::trace_back()$call, 2)[[1]]
    )
    class(tbp_object) <- "tbp"
    modeldata_temp[[f_name]] <- tbp_object
  }
  return(modeldata_temp)
}

#' Print to-be-provided object
#' @export
#' @keywords internal
print.tbp <- function(x) {
  req_data <- c()
  req_assumptions <- c()
  if (length(x$required_data) > 0) {
    req_data <- paste(
      "data:", paste(x$required_data, collapse = ",")
    )
  }
  if (length(x$required_assumptions) > 0) {
    req_assumptions <- paste(
      "assumptions:", paste(x$required_assumptions, collapse = ",")
    )
  }
  print(paste("Requires", paste(c(req_data, req_assumptions), collapse = " and ")))
}

#' @keywords internal
solve.tbp <- function(x, modeldata, data = list(), assumptions = list()) {
  provided <- x$check_provided(
    data = data, assumptions = assumptions,
    required_data = x$required_data, required_assumptions = x$required_assumptions
    )
  if (provided) {
    modeldata <- x$f_func(
      modeldata, data = data, assumptions = assumptions,
      required_data = x$required_data, required_assumptions = x$required_assumptions
      )
    modeldata[[x$f_name]] <- NULL
  }
  return(modeldata)
}

#' To be evaluated
#' @keywords internal
tbe <- function(r_expr, required = c(), calling_env = rlang::caller_env()) {
  if (modeldata_check(calling_env$modeldata, required, throw_error = FALSE)) {
    return(r_expr)
  } else {
    lazy_r <- lazyeval::lazy(r_expr)
    lazy_r$env <- calling_env

    tbe_object <- list()
    tbe_object <- list(
      lazy_r = lazy_r,
      required = required,
      calling_f = tail(rlang::trace_back()$call, 2)[[1]]
    )
    class(tbe_object) <- "tbe"
    return(tbe_object)

  }
}

#' Print to-be-evaluated object
#' @export
#' @keywords internal
print.tbe <- function(x) {
  print(paste("Waiting for information:", rlang::expr_text(x$lazy_r$expr)))
}

#' @keywords internal
solve.tbe <- function(x, modeldata, throw_error = TRUE) {
  evaluate <- modeldata_check(
    modeldata, x$required, calling_env = x$calling_f, throw_error = throw_error
  )
  if (evaluate) {
    return(lazyeval::lazy_eval(x$lazy_r, list(modeldata = modeldata)))
  } else {
    return(x)
  }
}

#' Lazy-execute a function once all necessary inputs are available
#'
#' @param f_name Name of the "functions" that will be executed once all
#'   necessary inputs are available.
#' @param f_expr Expression with arbitrary R code in which attributes are
#'   assigned to `modeldata.`
#' @param required A character vector with the names of required inputs.
#' @param modeldata Modeldata object in the calling function.
#' @param calling_env Calling environment, should be `rlang::caller_env()`.
#'
#' @details The difference between `tbp` (to be provided) and `tbc` (to be
#'   computed) is that `tbp` is for promises based on data and assumptions, and
#'   `tbc` is for promises based on internal modeldata variables (typically in
#'   `.metainfo`). This means that there can also be cascades: once a `tbp` is
#'   resolved, this may lead to an update of some metainformation, such that a
#'   `tbc` can be resolved in the next update.
#' @keywords internal
tbc <- function(f_name, f_expr, required = c(), required_values = NULL,
                advice = NULL,
                modeldata = NULL, calling_env = rlang::caller_env()) {
  if (is.null(modeldata)) {
    modeldata_temp <- calling_env$modeldata
  } else {
    modeldata_temp <- modeldata
  }

  f_lazy <- lazyeval::lazy(f_expr)
  f_lazy$env <- copy_env(f_lazy$env, exclude = "modeldata")

  f_func <- function(modeldata) {
    f_lazy$env$modeldata <- modeldata
    lazyeval::lazy_eval(f_lazy)
    return(f_lazy$env$modeldata)
  }

  environment(f_func) <- copy_env(calling_env, exclude = "modeldata")
  environment(f_func)$f_lazy <- f_lazy

  computable <- modeldata_check(
    modeldata_temp, required = required, required_values = required_values,
    advice = advice, throw_error = FALSE
  )

  if (computable) {
    modeldata_temp <- f_func(modeldata_temp)
  } else {
    tbc_object <- list(
      f_name = f_name,
      f_func = f_func,
      required = required,
      required_values = required_values,
      advice = advice,
      calling_f <- tail(rlang::trace_back()$call, 2)[[1]]
      )
    class(tbc_object) <- "tbc"
    modeldata_temp[[f_name]] <- tbc_object
  }
  return(modeldata_temp)
}

#' Print to-be-computed object
#' @export
#' @keywords internal
print.tbc <- function(x) {
  print(paste(
    "Waiting for the following information to compute values:",
    paste(x$required, collapse = ", ")
  ))
}

#' @keywords internal
solve.tbc <- function(x, modeldata, throw_error = TRUE) {
  computable <- modeldata_check(
    modeldata, x$required, required_values = x$required_values,
    advice = x$advice, calling_env = x$calling_f, throw_error = throw_error
  )
  if (computable) {
    modeldata <- x$f_func(modeldata)
    modeldata[[x$f_name]] <- NULL
  }
  return(modeldata)
}
