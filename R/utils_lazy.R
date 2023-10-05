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
  print(paste("Requires assumption:", lazyeval::expr_text(x$lazy_r$expr)))
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
  print(paste("Waiting for information:", lazyeval::expr_text(x$lazy_r$expr)))
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
