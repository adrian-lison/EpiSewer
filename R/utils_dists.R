# Exponential ----
## Discretized ----
get_discrete_exponential <- function(exponential_mean, maxX = NULL) {
  lambda <- 1 / exponential_mean

  # determine maxX
  if (is.null(maxX)) {
    maxX <- which(sapply(1:100, function(maxX) {
      (1 - pexp(maxX, rate = lambda)) < 0.005
    }))[1]
    if (is.na(maxX)) {
      maxX <- 100
      cli::cli_inform(c(
        "!" = paste0("Maximum length of distribution was set to 100. ",
                     "The last bin covers ",
                     100 * round(
                       (1 - pexp(maxX, rate = lambda)),
                       2
                     ),
                     "% of the probability mass."
        )
      ))
    }
  }

  # compute discrete distribution
  longest <- (1 - extraDistr::pdgamma(maxX, shape = 1, rate = lambda))
  probs <- c(
    # all except longest (discrete)
    extraDistr::ddgamma(0:(maxX - 1), shape = 1, rate = lambda),
    longest
  )
  return(probs)
}

# Gamma ----

## Parameterization ----
#' Get shape of a Gamma distribution given its mean and sd
#' @keywords internal
get_gamma_shape_alternative <- function(gamma_mean, gamma_sd) {
  gamma_shape <- (gamma_mean / gamma_sd)^2
  return(gamma_shape)
}

#' Get rate of a Gamma distribution given its mean and sd
#' @keywords internal
get_gamma_rate_alternative <- function(gamma_mean, gamma_sd) {
  gamma_rate <- gamma_mean / (gamma_sd^2)
  return(gamma_rate)
}

#' Get scale of a Gamma distribution given its mean and sd
#' @keywords internal
get_gamma_scale_alternative <- function(gamma_mean, gamma_sd) {
  return(1 / get_gamma_rate_alternative(gamma_mean, gamma_sd))
}

#' Get mean of a Gamma distribution given its shape and scale
#' @keywords internal
get_gamma_mean_alternative <- function(gamma_shape, gamma_scale) {
  gamma_mean <- gamma_shape * gamma_scale
  return(gamma_mean)
}

#' Get sd of a Gamma distribution given its shape and scale
#' @keywords internal
get_gamma_sd_alternative <- function(gamma_shape, gamma_scale) {
  gamma_sd <- sqrt(gamma_shape) * gamma_scale
  return(gamma_sd)
}

## Discretized ----

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
#' @param gamma_cv Alternative parameterization: Coefficient of variation of the
#'   Gamma distribution.
#' @param maxX Right truncation point. All probability mass beyond `maxX` will
#'   be assigned to `maxX`. If `NULL` (default), this is automatically chosen
#'   such that the last bin has less than 0.5% of the probability mass.
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
                               gamma_cv,
                               maxX = NULL,
                               include_zero = TRUE,
                               print_params = FALSE) {
  if (missing(gamma_shape)) {
    if (missing(gamma_mean) || (missing(gamma_sd) && missing(gamma_cv))) {
      stop("No valid combination of parameters supplied", call. = FALSE)
    }
    if (missing(gamma_sd)) {
      gamma_sd <- gamma_mean * gamma_cv
    }
    gamma_shape <- get_gamma_shape_alternative(gamma_mean, gamma_sd)
  }
  if (missing(gamma_rate)) {
    if (missing(gamma_scale)) {
      if (missing(gamma_mean) || (missing(gamma_sd) && missing(gamma_cv))) {
        stop("No valid combination of parameters supplied", call. = FALSE)
      }
      if (missing(gamma_sd)) {
        gamma_sd <- gamma_mean * gamma_cv
      }
      gamma_rate <- get_gamma_rate_alternative(gamma_mean, gamma_sd)
    } else {
      gamma_rate <- 1 / gamma_scale
    }
  }

  # shortest period (combines periods 0 and 1)
  shortest <- pgamma(2, shape = gamma_shape, rate = gamma_rate)
  # longest period (combines all periods >= maxX)
  if (is.null(maxX)) {
    maxX <- which(sapply(1:100, function(maxX) {
      (1 - pgamma(maxX, shape = gamma_shape, rate = gamma_rate)) < 0.005
    }))[1]
    if (is.na(maxX)) {
      maxX <- 100
      cli::cli_inform(c(
        "!" = paste0("Maximum length of distribution was set to 100. ",
                    "The last bin covers ",
                    100 * round(
                      (1 - pgamma(maxX, shape = gamma_shape, rate = gamma_rate)),
                      2
                      ),
                    "% of the probability mass."
        )
      ))
    }
  }
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
#'   be assigned to `maxX`.  If `NULL` (default), this is automatically chosen
#'   such that the last bin has approximately less than 0.5% of the probability
#'   mass.
#'
#' @return A numeric vector representing the PMF of the discretized shifted
#'   Gamma distribution.
#' @export
get_discrete_gamma_shifted <- function(
    gamma_mean, gamma_sd, maxX = NULL) {
  maxX <- 1 + length(get_discrete_gamma(
    gamma_mean = gamma_mean,
    gamma_sd = gamma_sd,
    maxX = maxX
  ))
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

# Log-Normal ----

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1 # error-function
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2) # inverse error function

## Parameterization ----

#' Compute mu parameter of Log-Normal from other quantities.
#'
#' @param unit_q5 5% quantile of distribution.
#' @param unit_q95 95% quantile of distribution.
#'
#' @details Currently, only conversion from quantiles is supported, but other
#' alternatives may be added.
#'
#' @return Mu parameter of Log-Normal distribution.
#' @keywords internal
get_lognormal_mu_alternative <- function(unit_mean = NULL, unit_sd = NULL, unit_q5 = NULL, unit_q95 = NULL) {
  if (!is.null(unit_mean) && !is.null(unit_sd)) {
    sigma2 <- log((unit_sd / unit_mean)^2 + 1)
    mu <- log(unit_mean) - sigma2 / 2
  } else if (!is.null(unit_q5) && !is.null(unit_q95)) {
    erfq5 <- erfinv(2*0.05-1)
    erfq95 <- erfinv(2*0.95-1)
    if (unit_q5 > unit_q95) {
      cli::cli_abort(paste(
        "Lower quantile `unit_q5` must not be larger",
        "than upper quantile `unit_q95`."
      ))
    }
    mu = (log(unit_q5)/erfq5 - log(unit_q95)/erfq95) / (1/erfq5 - 1/erfq95)
  } else {
    cli::cli_abort(paste(
      "Either `unit_mean` and `unit_sd` or `unit_q5` and `unit_q95` must be supplied."
    ))
  }
  return(mu)
}

#' Compute sigma parameter of Log-Normal from other quantities.
#'
#' @param mu Mu parameter of Log-Normal distribution.
#' @inheritParams get_lognormal_mu_alternative
#'
#' @details Currently, only conversion from mu + a quantile is supported, but
#'   other alternatives may be added.
#'
#' @return Sigma parameter of Log-Normal distribution.
#' @keywords internal
get_lognormal_sigma_alternative <- function(mu, unit_mean = NULL, unit_sd = NULL,
                                            unit_q5 = NULL, unit_q95 = NULL) {
  if (!is.null(unit_mean) && !is.null(unit_sd)) {
    sigma <- sqrt(log((unit_sd / unit_mean)^2 + 1))
  } else {
    if (is.null(unit_q5) && is.null(unit_q95)) {
      cli::cli_abort(
        "Either `unit_q5` or `unit_q95` must be supplied together with mu."
        )
    }
    if (!is.null(unit_q5) && !is.null(unit_q95)) {
      cli::cli_warn(paste(
        "Both `unit_q5` and `unit_q95` were supplied together with mu,",
        "using only `unit_q95` to compute sigma."
      ))
    }
    if (!is.null(unit_q95)) {
      if (unit_q95 < exp(mu)) {
        cli::cli_abort(
          "Upper quantile `unit_q95` must not be less than exp(mu)."
        )
      }
      sigma = (log(unit_q95) - mu)/(sqrt(2) * erfinv(2*0.95-1))
    } else if (!is.null(unit_q5)) {
      if (unit_q5 > exp(mu)) {
        cli::cli_abort(
          "Lower quantile `unit_q5` must not be greater than exp(mu)."
        )
      }
      sigma = (log(unit_q5) - mu)/(sqrt(2) * erfinv(2*0.05-1))
    }
  }
  return(sigma)
}

## Discretized ----

#' Get PMF of a discretized lognormal distribution
#'
#' This function accepts both log-scale and unit-scale parameters to specify a
#' discretized lognormal distribution.
#'
#' @param meanlog Log scale mean (location of lognormal).
#' @param sdlog Log scale standard deviation (scale of lognormal).
#' @param unit_mean Alternative parameterization: unit/natural scale mean.
#' @param unit_sd Alternative parameterization: unit/natural scale sd.
#' @param unit_cv Alternative parameterization: unit/natural scale coefficient
#'   of variation.
#' @param maxX Right truncation point. All probability mass beyond `maxX` will
#'   be assigned to `maxX`. If `NULL` (default), this is automatically chosen
#'   such that the last bin has less than 0.5% of the probability mass.
#' @param include_zero Should the distribution explicitly cover X=0, or should
#'   X=1 include the probability mass for X=0 too?
#' @param print_params Should the log-level parameters be printed?
#'
#' @return A numeric vector representing the PMF of the discretized lognormal.
#' @export
get_discrete_lognormal <- function(
    meanlog, sdlog, unit_mean = NULL, unit_sd = NULL, unit_cv = NULL, maxX = NULL, include_zero = TRUE,
    print_params = FALSE) {
  if (!is.null(unit_mean) && (!is.null(unit_sd) || !is.null(unit_cv))) {
    if (is.null(unit_sd)) {
      unit_sd <- unit_mean * unit_cv
    }
    sigma2 <- log((unit_sd / unit_mean)^2 + 1)
    meanlog <- log(unit_mean) - sigma2 / 2
    sdlog <- sqrt(sigma2)
  }
  # shortest period (combines periods 0 and 1)
  shortest <- plnorm(2, meanlog = meanlog, sdlog = sdlog)
  # longest period (combines all periods >= maxX)
  if (is.null(maxX)) {
    maxX <- which(sapply(1:100, function(maxX) {
      (1 - plnorm(maxX, meanlog = meanlog, sdlog = sdlog)) < 0.005
    }))[1]
    if (is.na(maxX)) {
      maxX <- 100
      cli::cli_inform(c(
        "!" = paste0("Maximum length of distribution was set to 100. ",
                    "The last bin covers ",
                    100 * round(
                      (1 - plnorm(maxX, meanlog = meanlog, sdlog = sdlog)),
                      2
                    ),
                    "% of the probability mass."
                    )
        ))
    }
  }
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

# Truncated normal ----

#' Calculate the mean of a truncated normal distribution (truncated below zero)
#' with mean mu and standard deviation sigma.
#'
#' @param mu Mean of the untruncated normal distribution
#' @param sigma Standard deviation of the untruncated normal distribution
#'
#' @return The mean of the truncated normal distribution
trunc_normal_mean <- function(mu, sigma) {
  if (sigma == 0) {
    return(mu)
  } else {
    alpha <- -mu / sigma
    phi_alpha <- exp(dnorm(alpha, log = TRUE))  # std_normal_lpdf
    Phi_alpha <- pnorm(alpha)                   # Phi_approx
    return(mu + sigma * phi_alpha / (1 - Phi_alpha))
  }
}


# Beta ----
get_beta_alpha_alternative <- function(beta_mean, beta_sd) {
  alpha <- beta_mean * (beta_mean * (1 - beta_mean) / beta_sd^2 - 1)
  return(alpha)
}

get_beta_beta_alternative <- function(beta_mean, beta_sd) {
  beta <- (1 - beta_mean) * (beta_mean * (1 - beta_mean) / beta_sd^2 - 1)
  return(beta)
}

check_beta_alternative <- function(beta_mean, beta_sd) {
  if (!(beta_mean %in% c(0,1) && beta_sd == 0)) {
    alpha <- get_beta_alpha_alternative(beta_mean, beta_sd)
    beta <- get_beta_beta_alternative(beta_mean, beta_sd)
    if (alpha <= 0 || beta <= 0 || beta_sd < 0) {
      cli::cli_abort(paste(
        "Invalid beta distribution parameters supplied.",
        "The mean must be between 0 and 1, and the standard deviation",
        "must be positive but not too large."
      ))
    }
  }
}

#' Exponential-Gamma distribution
#'
#' The Exponential-Gamma (EG) distribution is also known as the Lomax
#' distribution and has a characteristic long tail.
#'
#' @param p Vector of probabilities.
#' @param shape Shape of the Exponential-Gamma distribution.
#' @param scale Scale of the Exponential-Gamma distribution.
qexpgamma <- function(p, shape, scale) {
  names(p) <- paste0(100*p,"%")
  round(sapply(p, extraDistr::qlomax, kappa = shape, lambda = 1/scale),3)
}

# Distribution validation ----

check_dist <- function(dist, name = "probability distribution") {
  if (!is.numeric(dist)) {
    cli::cli_abort(paste("Supplied", name, "is not a numeric vector."))
  }
  if (any(dist < 0)) {
    cli::cli_abort(paste(
      "Supplied", name, "has negative entries.",
      "All probabilities must be positive."
    ))
  }
  if (sum(dist) != 1) {
    cli::cli_warn(paste(
      "Supplied", name, "does not sum to 1.",
      "EpiSewer will normalize the probabilities such that they sum to 1.\n"
    ))
    dist <- dist / sum(dist)
  }
  return(dist)
}


#' Get the mean of a discrete distribution
#'
#' @param dist Discrete distribution represented as numeric vector
#' @param include_zero If `TRUE` (default), the vector index starts at zero.
#'   Otherwise, it starts at 1.
#'
#' @return Mean of the discrete distribution
#' @keywords internal
dist_mean <- function(dist, include_zero = TRUE) {
  if (include_zero) {
    mean <- sum((0:(length(dist) - 1)) * dist)
  } else {
    mean <- sum((1:(length(dist))) * dist)
  }
  return(mean)
}
