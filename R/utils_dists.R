## ---------------------------------------------------------------
##                Discretized delay distributions               -
## ---------------------------------------------------------------

#' Get shape of a Gamma distribution given its mean and sd
get_gamma_shape_alternative <- function(gamma_mean, gamma_sd) {
  gamma_shape <- (gamma_mean / gamma_sd)^2
  return(gamma_shape)
}

#' Get rate of a Gamma distribution given its mean and sd
get_gamma_rate_alternative <- function(gamma_mean, gamma_sd) {
  gamma_rate <- gamma_mean / (gamma_sd^2)
  return(gamma_rate)
}

#' Get scale of a Gamma distribution given its mean and sd
get_gamma_scale_alternative <- function(gamma_mean, gamma_sd) {
  return(1 / get_gamma_rate_alternative(gamma_mean, gamma_sd))
}

#' Get PMF of a discretized Gamma distribution
#'
#' @description This function accepts different parameterizations to specify a
#'   discretized Gamma distribution.
#'
#' @param gamma_shape Shape parameter of the Gamma distribution.
#' @param gamma_rate Rate parameter of the Gamma distribution.
#' @param gamma_scale Scale parameter of the Gamma distribution. Can be
#'   specified instead of the rate. Only has an effect if the rate is not
#'   specified.
#' @param gamma_mean Alternative parameterization: Mean of the Gamma
#'   distribution.
#' @param gamma_sd Alternative parameterization: Standard deviation of the Gamma
#'   distribution.
#' @param maxX Right truncation point. All probability mass beyond `maxX` will
#'   be assigned to `maxX`.
#' @param include_zero Should the distribution explicitly cover X=0, or should
#'   X=1 include the probability mass for X=0 too?
#' @param print_params Should the shape and rate parameters be printed?
#'
#' @return A numeric vector representing the PMF of the discretized Gamma
#'   distribution.
#' @export
get_discrete_gamma <- function(gamma_shape,
                               gamma_rate,
                               gamma_scale,
                               gamma_mean,
                               gamma_sd,
                               maxX,
                               include_zero = TRUE,
                               print_params = FALSE) {
  if (missing(gamma_shape)) {
    if (missing(gamma_mean) || missing(gamma_sd)) {
      stop("No valid combination of parameters supplied", call. = FALSE)
    }
    gamma_shape <- get_gamma_shape_alternative(gamma_mean, gamma_sd)
  }
  if (missing(gamma_rate)) {
    if (missing(gamma_scale)) {
      if (missing(gamma_mean) || missing(gamma_sd)) {
        stop("No valid combination of parameters supplied", call. = FALSE)
      }
      gamma_rate <- get_gamma_rate_alternative(gamma_mean, gamma_sd)
    } else {
      gamma_rate <- 1 / gamma_scale
    }
  }

  # shortest period (combines periods 0 and 1)
  shortest <- pgamma(2, shape = gamma_shape, rate = gamma_rate)
  # longest period (combines all periods >= maxX)
  longest <- (1 - pgamma(maxX, shape = gamma_shape, rate = gamma_rate))

  if (include_zero) {
    probs <- c(
      # all except longest (discrete)
      extraDistr::ddgamma(0:(maxX - 1), shape = gamma_shape, rate = gamma_rate),
      longest
    )
  } else {
    probs <- c(
      shortest,
      # all except shortest and longest (discrete)
      extraDistr::ddgamma(2:(maxX - 1), shape = gamma_shape, rate = gamma_rate),
      longest
    )
  }

  if (print_params) {
    print(paste("Shape =", gamma_shape, "| Rate =", gamma_rate))
  }

  return(probs)
}

#' Get PMF of a discretized shifted Gamma distribution
#'
#' @description This function specifies a discretized shifted Gamma distribution
#'   (shifted such that the minimum is at 1) using a mean and standard deviation
#'   parameterization.
#'
#' @description The shift makes the distribution attractive for modeling
#'   generation time distributions, which are often assumed to have zero
#'   probability mass for a generation time of zero (as this is incompatible
#'   with the renewal equation).
#'
#' @details This code was adapted from EpiEstim, credit to Anne Cori
#'   (a.cori@imperial.ac.uk).
#'
#' @param gamma_mean Mean of the shifted Gamma distribution.
#' @param gamma_sd Standard deviation of the shifted Gamma distribution.
#' @param maxX Right truncation point. All probability mass beyond `maxX` will
#'   be assigned to `maxX`.
#'
#' @return A numeric vector representing the PMF of the discretized shifted
#'   Gamma distribution.
#' @export
get_discrete_gamma_shifted <- function(
    gamma_mean = 4.8, gamma_sd = 2.3, maxX = 10) {
  k <- 1:maxX
  if (gamma_sd < 0) {
    stop("gamma_sd must be >=0.")
  }
  if (gamma_mean <= 1) {
    stop("gamma_mean must be >1")
  }
  if (maxX < 1) {
    stop("Maximum period must be >=1.")
  }
  a <- ((gamma_mean - 1) / gamma_sd)^2
  b <- gamma_sd^2 / (gamma_mean - 1)
  cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, scale = b)
  res <- k * cdf_gamma(k, a, b) + (k - 2) *
    cdf_gamma(k - 2, a, b) - 2 * (k - 1) * cdf_gamma(k - 1, a, b)
  res <- res + a * b * (2 * cdf_gamma(k - 1, a + 1, b) -
    cdf_gamma(k - 2, a + 1, b) - cdf_gamma(k, a + 1, b))
  res <- sapply(res, function(e) max(0, e))
  res <- res / sum(res)
  return(res)
}

#' Get PMF of a discretized lognormal distribution
#'
#' This function accepts both log-scale and unit-scale parameters to specify a
#' discretized lognormal distribution.
#'
#' @param meanlog Log scale mean (location of lognormal).
#' @param sdlog Log scale standard deviation (scale of lognormal).
#' @param unit_mean Alternative parameterization: unit scale mean.
#' @param unit_sd Alternative parameterization: unit scale sd.
#' @param maxX Right truncation point. All probability mass beyond `maxX` will
#'   be assigned to `maxX`.
#' @param include_zero Should the distribution explicitly cover X=0, or should
#' X=1 include the probability mass for X=0 too?
#' @param print_params Should the log-level parameters be printed?
#'
#' @return A numeric vector representing the PMF of the discretized lognormal.
#' @export
get_discrete_lognormal <- function(
    meanlog, sdlog, unit_mean = NULL, unit_sd = NULL, maxX, include_zero = TRUE,
    print_params = FALSE) {
  if (!is.null(unit_mean) && !is.null(unit_sd)) {
    sigma2 <- log((unit_sd / unit_mean)^2 + 1)
    meanlog <- log(unit_mean) - sigma2 / 2
    sdlog <- sqrt(sigma2)
  }
  # shortest period (combines periods 0 and 1)
  shortest <- plnorm(2, meanlog = meanlog, sdlog = sdlog)
  # longest period (combines all periods >= maxX)
  longest <- (1 - plnorm(maxX, meanlog = meanlog, sdlog = sdlog))

  if (include_zero) {
    # all except longest (discrete)
    probs <- c(
      sapply(0:(maxX - 1), function(x) {
        plnorm(
          x + 1,
          meanlog = meanlog, sdlog = sdlog
        ) -
          plnorm(x, meanlog = meanlog, sdlog = sdlog)
      }),
      longest
    )
  } else {
    probs <- c(
      shortest,
      # all except shortest and longest (discrete)
      sapply(2:(maxX - 1), function(x) {
        plnorm(
          x + 1,
          meanlog = meanlog, sdlog = sdlog
        ) -
          plnorm(x, meanlog = meanlog, sdlog = sdlog)
      }),
      longest
    )
  }

  if (print_params) {
    print(paste("meanlog =", meanlog, "| sdlog =", sdlog))
  }

  return(probs)
}

check_dist <- function(dist, name = "probability distribution") {
  if (!is.numeric(dist)) {
    rlang::abort(paste("Supplied", name, "is not a numeric vector."))
  }
  if (any(dist < 0)) {
    rlang::abort(paste(
      "Supplied", name, "has negative entries.",
      "All probabilities must be positive."
    ))
  }
  if (sum(dist) != 1) {
    rlang::warn(paste(
      "Supplied", name, "does not sum to 1.",
      "EpiSewer will normalize the probabilities such that they sum to 1.\n"
    ))
    dist <- dist / sum(dist)
  }
  return(dist)
}
