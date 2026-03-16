suppress_warnings <- function(.expr, .f, ...) {
  eval.parent(substitute(
    withCallingHandlers(.expr, warning = function(w) {
      cm <- conditionMessage(w)
      cond <- any(sapply(.f_list, function(.f) {
        if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
      }
      ))
      if (cond) {
        invokeRestart("muffleWarning")
      }
    })
  ))
}

suppress_messages <- function(.expr, .f_list, ...) {
  eval.parent(substitute(
    withCallingHandlers(
      .expr,
      message = function(w) {
        cm <- conditionMessage(w)
        cond <- any(sapply(.f_list, function(.f) {
          if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
        }
        ))
        if (cond) {
          invokeRestart("muffleMessage")
        }
      }
    )
  ))
}

suppress_messages_warnings <- function(.expr, .f_list, ...) {
  eval.parent(substitute(
    withCallingHandlers(
      .expr,
      message = function(w) {
        cm <- conditionMessage(w)
        cond <- any(sapply(.f_list, function(.f) {
          if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
        }
        ))
        if (cond) {
          invokeRestart("muffleMessage")
        }
      },
      warning = function(w) {
        cm <- conditionMessage(w)
        cond <- any(sapply(.f_list, function(.f) {
          if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
        }
        ))
        if (cond) {
          invokeRestart("muffleWarning")
        }
      }
    )
  ))
}

#' Semantic sugar to wrap an abort statement in a function
#'
#' This can be used to provide abort statements as a direct argument to
#' tryCatch without wrapping them in a function call
#' @keywords internal
abort_f <- function(...) {
  return(function(err) {
    cli::cli_abort(
      call = NULL, .trace_bottom = sys.frames()[[sys.nframe() - 4]], ...
    )
  })
}

#' Collect warnings arising during computation of a result
#'
#' @param expr Expression with the computation to run
#'
#' @return A `list`, with two elements: value and warnings
#' @keywords internal
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}
