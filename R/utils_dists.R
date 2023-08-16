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

#' Get PMF of a discretized Gamma distribution.
#'
#' This function accepts different parameterizations to specify the Gamma
#' distribution
#'
#' @param gamma_shape Shape parameter of the Gamma distribution
#' @param gamma_rate Rate parameter of the Gamma distribution.
#' @param gamma_scale Scale parameter of the Gamma distribution. Can be
#'   specified instead of the rate. Only has an effect if the rate is not
#'   specified.
#' @param gamma_mean Alternative parameterization: Mean of the Gamma
#' @param gamma_sd Alternative parameterization: Standard deviation of the Gamma
#' @param maxX Right truncation point
#' @param include_zero Should the distribution explicitly cover X=0, or should
#'   X=1 include the probability mass for X=0 too?
#' @param print_params Should the shape and rate parameters be printed?
#'
#' @return PMF of the discretized Gamma distribution
#' @export
get_discrete_gamma <- function(gamma_shape,
                               gamma_rate,
                               gamma_scale,
                               gamma_mean,
                               gamma_sd,
                               maxX,
                               include_zero = T,
                               print_params = F) {
  if (missing(gamma_shape)) {
    if (missing(gamma_mean) || missing(gamma_sd)) {
      stop("No valid combination of parameters supplied", call. = F)
    }
    gamma_shape <- get_gamma_shape_alternative(gamma_mean, gamma_sd)
  }
  if (missing(gamma_rate)) {
    if (missing(gamma_scale)) {
      if (missing(gamma_mean) || missing(gamma_sd)) {
        stop("No valid combination of parameters supplied", call. = F)
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

#' Title
#'
#' @description # This ocode is adapted from EpiEstim, credit to Anne Cori
#'   a.cori@imperial.ac.uk
#'
#' @param gamma_mean
#' @param gamma_sd
#' @param maxX
#'
#' @return
#' @export
#'
#' @examples
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
  return(res)
}

#' Get PMF of a discretized lognormal distribution.
#'
#' This function accepts both log-scale and unit-scale parameters to specify the
#' lognormal distribution
#'
#' @param meanlog Mean of log
#' @param sdlog Standard deviation of log
#' @param unit_mean Alternative parameterization: unit scale mean
#' @param unit_sd Alternative parameterization: unit scale sd
#' @param maxX Right truncation point
#' @param include_zero Should the distribution explicitly cover X=0, or should
#' X=1 include the probability mass for X=0 too?
#' @param print_params Should the log-level parameters be printed?
#'
#' @return PMF of the discretized lognormal
#' @export
get_discrete_lognormal <- function(
    meanlog, sdlog, unit_mean = NULL, unit_sd = NULL, maxX, include_zero = T,
    print_params = F) {
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
