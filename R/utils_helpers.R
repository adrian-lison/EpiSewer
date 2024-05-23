## ----------------------------------------------------------------
##                        Helper functions                       -
## ----------------------------------------------------------------

#' Clip values of a vector between a lower and an upper bound
#'
#' @param vec A vector with values to be clipped
#' @param LB Lower bound
#' @param UB Upper bound
#'
#' @return A vector with clipped/fenced values
fence <- function(vec, LB = -Inf, UB = Inf) {
  pmax(LB, pmin(vec, UB))
}

#' Compare entries of two vectors with potentially missing values
#'
#' @param v1 First vector
#' @param v2 Second vector, should be of same length as first one
#'
#' @return A boolean vector with element-wise comparisons. Comparisons including
#' NAs are coded as FALSE.
compareNA <- function(v1, v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

#' Compute a weighted median of values
#'
#' @param x A vector with values to be aggregated
#' @param w Weights for each entry
#'
#' @return Weighted median of the values in x
weighted.median <- function(x, w) {
  w <- w[order(x)]
  x <- x[order(x)]

  prob <- cumsum(w) / sum(w)
  ps <- which.min(abs(prob - .5))
  return(x[ps])
}

#' Vectorized version of `rep` with `each` argument, `each` can be a vector
rep_each_v <- function(x, each) {
  unlist(mapply(function(i, e) rep(i, each = e), i = x, e = each))
}


#' Places knots for fitting B-splines to a time series
#'
#' @param ts_length Length of the time series
#' @param knot_distance Normal distance between knots
#' @param partial_window Window with only partial data towards the present in
#'   which the knot distances will be shorter. This is to avoid erroneous
#'   extrapolation of splines towards the present.
#'
#' @return A vector with knot positions.
place_knots <- function(ts_length, knot_distance, partial_window = 30) {
  if (!knot_distance>1) {
    cli::cli_abort("Knot distance must be larger than one.")
  }
  if (!ts_length>knot_distance) {
    cli::cli_abort("The length of the time series must be larger than the knot distance.")
  }
  if (!partial_window>0) {
    cli::cli_abort("The `partial_window` must be larger than zero.")
  }
  # define knot distances close to present, i.e. in window with partial data
  last_dists <- 1+2^seq(1, ceiling(log2(knot_distance)))
  last_dists <- last_dists[last_dists<=knot_distance]
  last_dists <- c(
    rep(3, max(ceiling((partial_window-sum(last_dists))/3),1)),
    last_dists
    )
  last_dists <- last_dists[1:which(cumsum(last_dists)>=partial_window)[1]]
  last_dists <- last_dists[cumsum(last_dists)<ts_length+1]
  # define knot distances for remaining window (full data)
  remaining_knots <- (ts_length+1-sum(last_dists)) %/% knot_distance
  if (remaining_knots > 0) {
    all_dists <- c(last_dists, rep(knot_distance, remaining_knots))
  } else {
    all_dists <- last_dists
  }
  # get knot positions
  int_knots <- rev(ts_length+1-cumsum(all_dists))
  int_knots <- c(min(int_knots) - diff(int_knots)[1], int_knots) # add one internal knot before t=0
  bound_knots <- c(min(int_knots) - diff(int_knots)[1], ts_length + 1)
  return(list(interior = int_knots, boundary = bound_knots))
}


#' Obtain regression vector to estimate a linear trend
#'
#' @description This function provides the solution vector from a (weighted)
#'   linear regression with one independent variable x, where x may represent
#'   points in a time series.
#'
#' @param x A `vector` with observations of the single independent variable (for
#'   example points in time).
#' @param weights Optional, a `vector` with weights for each observation.
#'
#' @return A `vector` which contains the relevant row of the regression solution
#'   matrix that represents the trend coefficient. By multiplying this vector
#'   with the observed time series values y, you directly get the (weighted
#'   least squares) estimate of the trend.
#'
#' @examples
#' x = 1:6
#' y = x*4 + rnorm(6, 0, 0.1)
#' w <- c(0.1, 0.1, 0.1, 0.2, 0.2, 0.3) # weights
#' trend_reg <- get_regression_linear_trend(x, weights = w)
#' as.vector(trend_reg %*% y) # this gives you the trend estimate
#' summary(lm(y ~ x, weights = w))$coefficients["x","Estimate"] # should be the same
get_regression_linear_trend <- function(x, weights = rep(1/length(x), length(x))) {
  if (length(weights) != length(x)) {
    cli::cli_abort(
      "Provided weights do not match number of observations.", .internal = TRUE
      )
  }
  X <- matrix(c(rep(1,length(x)), x), byrow = F, ncol = 2)
  W <- diag(weights)
  trend_reg <- (solve(t(X) %*% W %*% X) %*% t(X) %*% W)[2,]
  return(trend_reg)
}

logistic_find_k <- function(a, c) {
  eq <- function(k) {
    (a * c * k * exp(k)) / ((a + exp(k))^2) - 1
  }
  result <- uniroot(eq, interval = c(0.9, 2))
  return(result$root)
}

logistic_find_a_k <- function(c) {
  k_value <- 1
  k_value_old <- 0
  i <- 1
  while (i<1000 && abs(k_value_old-k_value)>0.0001) {
    a_value <- (c-1)/exp(-k_value)
    k_value_old <- k_value
    k_value <- logistic_find_k(a_value, c)
    i <- i + 1
  }
  return(list(a = a_value, k = k_value))
}

logistic <- function(x, c, a, k) {
  return(c / (1+a*exp(-k*x)))
}

logistic_deriv <- function(x, c, a, k) {
  return(a*c*k*exp(k*x)/((a+exp(k*x))^2))
}

check_list_nested <- function(list_to_check, flat_var, var_vals = NULL) {
  if (is.null(var_vals)) {
    var_vals <- rep(NA, length(flat_var))
  }
  vars_levels <- stringr::str_split(flat_var, "\\$")
  presence <- mapply(function(levels, value) {
    check_l <- list_to_check
    for (l in levels) {
      if (!(l %in% names(check_l))) {
        return(FALSE)
      } else {
        check_l <- check_l[[l]]
      }
    }
    if (any(c("tbe", "tbc", "tbp") %in% class(check_l))) {
      return(FALSE)
    } else if (!is.na(value)) {
      return(check_l == value)
    } else {
      return(TRUE)
    }
  }, levels = vars_levels, value = var_vals)
  return(presence)
}

default_list_nested <- function(list_to_validate, levels, default, i = 1) {
  if (i < length(levels)) {
    if (!(levels[i] %in% names(list_to_validate))) {
      list_to_validate[[levels[i]]] <- default_list_nested(
        list(), levels, default, i + 1
      )
    } else {
      list_to_validate[[levels[i]]] <- default_list_nested(
        list_to_validate[[levels[i]]], levels, default, i + 1
      )
    }
  } else {
    if (!(levels[i] %in% names(list_to_validate))) {
      list_to_validate[[levels[i]]] <- default
    }
  }
  return(list_to_validate)
}

list_except <- function(l, except) {
  sel <- names(l)[!(names(l) %in% except)]
  return(l[sel])
}

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
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}
