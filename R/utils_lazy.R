#' To be provided
tbp <- function(f_name, f_expr,
                required_data = c(), required_assumptions = c(),
                modeldata = NULL,
                calling_env = rlang::caller_env()) {

  if (is.null(modeldata)) {
    calling_modeldata <- calling_env$modeldata
  } else {
    calling_modeldata <- modeldata
  }

  f_lazy <- lazyeval::lazy(f_expr)
  env_merge <- function(f_env = list(), data = list(), assumptions = list()) {
    for (r in required_data) {
      if (!(r %in% names(f_env) && !is.null(f_env[[r]]))) {
        f_env[[r]] <- data[[r]]
      }
    }
    for (r in required_assumptions) {
      if (!(r %in% names(f_env) && !is.null(f_env[[r]]))) {
        f_env[[r]] <- assumptions[[r]]
      }
    }
    return(f_env)
  }
  f_check <- function(f_env = list(), data = list(), assumptions = list()) {
    f_env <- env_merge(f_env, data, assumptions)
    if (all(c(required_data, required_assumptions) %in% names(f_env))) {
      return (!any(sapply(
        as.list(f_env)[c(required_data, required_assumptions)],
        is.null
        )))
    } else {
      return (FALSE)
    }
  }
  f_func <- function(modeldata, data = list(), assumptions = list()) {
    f_lazy$env <- env_merge(f_lazy$env, data, assumptions)
    f_lazy$env$modeldata <- modeldata
    lazyeval::lazy_eval(f_lazy) # apply function to modeldata
    f_sewer_data <- list()
    f_sewer_data[[f_name]] <- as.list(f_lazy$env)[required_data]
    f_lazy$env$modeldata$.sewer_data <- c(
      f_lazy$env$modeldata$.sewer_data, f_sewer_data
      )
    f_sewer_assumptions <- list()
    f_sewer_assumptions[[f_name]] <- as.list(f_lazy$env)[required_assumptions]
    f_lazy$env$modeldata$.sewer_assumptions <- c(
      f_lazy$env$modeldata$.sewer_assumptions, f_sewer_assumptions
      )
    return(f_lazy$env$modeldata)
  }
  rm(modeldata, envir = calling_env)
  calling_env$f_name <- f_name
  calling_env$f_lazy <- f_lazy
  calling_env$env_merge <- env_merge
  calling_env$required_data = required_data
  calling_env$required_assumptions = required_assumptions
  environment(f_func) <- calling_env

  if (f_check(calling_env)) {
    calling_modeldata <- f_func(calling_modeldata)

  } else {
    tbp_o <- list(
      name = f_name,
      func = f_func,
      check = f_check,
      required_data = required_data,
      required_assumptions = required_assumptions,
      calling_f = tail(rlang::trace_back()$call, 2)[[1]])
    class(tbp_o) <- "tbp"
    calling_modeldata[[f_name]] <- tbp_o
    return(calling_modeldata)
  }
}

#' Print to-be-provided object
#' @export
print.tbp <- function(x) {
  req_data = c()
  req_assumptions = c()
  if (length(x$required_data) > 0) {
    req_data = paste(
      "data:", paste(x$required_data, collapse = ",")
      )
  }
  if (length(x$required_assumptions) > 0) {
    req_assumptions = paste(
      "assumptions:", paste(x$required_assumptions, collapse = ",")
      )
  }
  print(paste("Requires", paste(c(req_data, req_assumptions), collapse = " and ")))
}

solve.tbp <- function(x, modeldata, data = list(), assumptions = list()) {
  if (x$check(data = data, assumptions = assumptions)) {
    modeldata <- x$func(modeldata, data, assumptions)
    modeldata[[x$name]] <- NULL
  }
  return(modeldata)
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

#' Print to-be-evaluated object
#' @export
print.tbe <- function(x) {
  print(paste("Waiting for information:", rlang::expr_text(x$lazy_r$expr)))
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
tbc <- function(f_name, f_expr, required = c(), modeldata = NULL,
                calling_env = rlang::caller_env()) {
  if (is.null(modeldata)) {
    calling_modeldata <- calling_env$modeldata
  } else {
    calling_modeldata <- modeldata
  }

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

#' Print to-be-computed object
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
