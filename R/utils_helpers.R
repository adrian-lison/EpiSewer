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

check_list_nested <- function(list_to_check, flat_var) {
  vars_levels <- stringr::str_split(flat_var,"\\$")
  presence <- sapply(vars_levels, function(levels) {
    check_l <- list_to_check
    for (l in levels) {
      if (!(l %in% names(check_l))) {
        return(FALSE)
      } else {
        check_l <- check_l[[l]]
      }
    }
    if (any(c("tbe","tbef") %in% class(check_l))) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
return(presence)
}

default_list_nested <- function(list_to_validate, levels, default, i = 1) {
  if (i < length(levels)) {
    if (!(levels[i] %in% names(list_to_validate))) {
      list_to_validate[[levels[i]]] <- default_list_nested(list(), levels, default, i+1)
    } else {
      list_to_validate[[levels[i]]] <- default_list_nested(list_to_validate[[levels[i]]], levels, default, i+1)
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
