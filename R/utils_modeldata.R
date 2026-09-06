#' Template function
#'
#' @description This function does not actually do anything. It only serves as
#'   a template for the documentation of other functions using inheritance.
#'
#' @return A `model_component` object describing this part of the model
#'   specification. Pass it to the corresponding module function (e.g.
#'   [model_measurements()], [model_infections()]), which in turn is supplied
#'   to [EpiSewer()].
#' @keywords internal
template_model_helpers <- function() { }

#' Get modeling functions for a component (internal function)
#'
#' Get all available modeling functions for a given component. Argument defaults
#' are such that helpers are returned as a comma-separated list of roxygen
#' function links.
#'
#' @param component A character vector with the name of the component
#' @param collapse A character used to collapse helpers into one string. Default
#'   is NULL (do not collapse).
#' @param prefix A prefix to add to each function name
#' @param suffix A suffix to add to each function name
#'
#' @return A character vector with all available modeling functions for the
#'   component. If collapse is not NULL, the helper functions are collapsed into
#'   a single string.
#' @keywords internal
component_functions_ <- function(component, collapse = "\n",
                                 prefix = "- [", suffix = "()]") {
  if (!component %in% all_components()) {
    err_message <- c(
      paste(
        "No valid component provided.",
        "Must be one out of:"
      ),
      all_components()
    )
    names(err_message) <- c("!", rep("*", length(err_message) - 1))
    cli::cli_abort(err_message)
  }
  all_fs <- names(rlang::ns_env("EpiSewer"))
  arg_fs <- all_fs[stringr::str_detect(
    # this excludes functions ending with _
    all_fs, paste0("^", component, "_", ".*(?<!_)$")
  )]

  # sort values in arg_fs according to following scheme
  # 1. values ending with "_none" first
  # 2. values ending with "_observe" second
  # 3. values ending with "_assume" third
  # 4. values ending with "_estimate" fourth
  # 5. values ending with "_estimate" plus another suffix fifth
  arg_fs_none <- arg_fs[stringr::str_detect(arg_fs, "_none$")]
  arg_fs_observe <- arg_fs[stringr::str_detect(arg_fs, "_observe$")]
  arg_fs_assume <- arg_fs[stringr::str_detect(arg_fs, "_assume$")]
  arg_fs_estimate <- arg_fs[stringr::str_detect(arg_fs, "_estimate$")]
  arg_fs_estimate_other <- arg_fs[stringr::str_detect(arg_fs, "_estimate_")]
  arg_fs_estimate_other <- arg_fs_estimate_other[
    !arg_fs_estimate_other %in% c(
      arg_fs_none, arg_fs_observe, arg_fs_assume, arg_fs_estimate
    )
  ]
  arg_fs_sorted <- c(
    arg_fs_none, arg_fs_observe, arg_fs_assume, arg_fs_estimate, arg_fs_estimate_other
  )

  if (length(arg_fs_sorted) > 0) {
    helpers <- paste0(prefix, arg_fs_sorted, suffix)
  } else {
    helpers <- c()
  }
  if (!is.null(collapse)) {
    return(paste(helpers, collapse = collapse))
  } else {
    return(helpers)
  }
}

#' Get a list of modeling functions for a component
#'
#' @description This function returns a list of all available modeling functions
#'   for a given component.
#'
#' @inheritParams component_functions_
#'
#' @return A character vector with all available modeling functions for the
#'   component.
#' @export
component_functions <- function(component) {
  return(component_functions_(
    component,
    collapse = NULL, prefix = "", suffix = "()"
  ))
}

#' Define a prior in modeldata
#' @keywords internal
set_prior <- function(param, dist = "normal", ...) {
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


#' Define a normal prior in modeldata
#'
#' @param param Name of the parameter for which the prior is defined.
#' @param mu Mean.
#' @param sigma Standard deviation.
#' @param two_sigma Two times the standard deviation. Useful for defining priors
#'   via the two-sigma rule-of-thumb (approximately 95% of probability mass).
#' @param q5 Lower quantile (5%).
#' @param q95 Upper quantile (95%).
#'
#' @return Prior specification for modeldata.
#' @keywords internal
set_prior_normal <- function(param, mu = NULL, sigma = NULL, two_sigma = NULL, q5 = NULL, q95 = NULL) {
  if (!is.null(q5) && !is.null(q95)) {
    if (q5 > q95) {
      cli::cli_abort(
        "The lower quantile (q5) must be smaller than the upper quantile (q95)",
      )
    }
    mu <- (q5 + q95) / 2
    sigma <- (q95 - q5) / (2 * qnorm(0.95))
  } else {
    if (is.null(mu)) {
      cli::cli_abort(
        "mu must be supplied for the normal prior",
        .internal = TRUE
      )
    }
    if (is.null(sigma)) {
      if (is.null(two_sigma)) {
        cli::cli_abort(
          "Either sigma or two_sigma must be supplied for the normal prior",
          .internal = TRUE
        )
      } else {
        sigma <- two_sigma / 2
      }
    }
  }
  return(set_prior(param = param, dist = "normal", mu = mu, sigma = sigma))
}

#' Define a normal prior on the log scale in modeldata using median and factor
#'
#' @description This parameterization is useful to define a best guess for the
#'   value of a parameter (median), and a maximum factor by which we expect to
#'   deviate from that guess.
#'
#' @param param Name of the parameter for which the prior is defined.
#' @param unit_median Median on the unit/natural scale.
#' @param unit_q5 5% quantile of the distribution on the unit scale.
#' @param unit_q95 95% quantile of the distribution on the unit scale.
#' @param unit_factor By which factor do we expect the true parameter value to
#'   differ at most from our prior median? For example, `unit_factor = 2` would
#'   mean that we expect the true parameter to be at most twice our prior
#'   median, and at least half of our prior median. This uses the two-sigma
#'   rule-of-thumb. Must be specified together with `unit_median.`
#' @param paste_dist Additional information that should be pasted to the dist
#'   description.
#'
#' @details A normal prior on the log scale is effectively a log-normal prior on
#'   the unit/natural scale.
#'
#' @return Prior specification for modeldata.
#' @keywords internal
set_prior_normal_log <- function(param,
                                 unit_median = NULL,
                                 unit_q5 = NULL, unit_q95 = NULL,
                                 unit_factor = NULL, paste_dist = "") {
  if (!is.null(unit_median)) {
    mu = log(unit_median)
    if (!is.null(unit_q5) || !is.null(unit_q95)) {
      sigma = get_lognormal_sigma_alternative(
        mu = mu, unit_q5 = unit_q5, unit_q95 = unit_q95
        )
    } else if (!is.null(unit_factor)) {
      sigma = log(unit_factor)/2
    } else {
      cli::cli_abort(paste(
        "You need to specify one additional argument besides",
        "`unit_median` to define the prior."
        ))
    }
  } else {
    if (!is.null(unit_q5) && !is.null(unit_q95)) {
      mu = get_lognormal_mu_alternative(unit_q5 = unit_q5, unit_q95 = unit_q95)
      sigma = get_lognormal_sigma_alternative(mu = mu, unit_q95 = unit_q95)
    } else {
      cli::cli_abort(paste(
        "Please specify either `unit_median`, or both `unit_q5` and",
        "`unit_q95` to define the prior."
      ))
    }
  }
  prior = set_prior(
    param = param, dist = paste0("normal", paste_dist), mu = mu , sigma = sigma
    )
  return(prior)
}

#' Define a truncated normal prior in modeldata
#'
#' @description This function defines a truncated normal prior (truncation at
#'   zero) for a parameter in modeldata.
#'
#' @param param Name of the parameter for which the prior is defined.
#' @param mu Mean (not accounting for truncation).
#' @param sigma Standard deviation (not accounting for truncation).
#' @param two_sigma Two times the standard deviation (not accounting for
#'   truncation). Useful for defining priors via the two-sigma rule-of-thumb
#'   (approximately 95% of probability mass).
#' @param q5 Lower quantile (5%).
#' @param q95 Upper quantile (95%).
#'
#' @return Prior specification for modeldata.
#' @keywords internal
set_prior_trunc_normal <- function(param, mu = NULL, sigma = NULL, two_sigma = NULL, q5 = NULL, q95 = NULL) {
  if (!is.null(q5) && !is.null(q95)) {
    if (q5 > q95) {
      cli::cli_abort(
        "The lower quantile (q5) must be smaller than the upper quantile (q95)",
      )
    }
    params <- find_mu_sigma_truncnorm(q5_target = q5, q95_target = q95)
    mu <- params$mu
    sigma <- params$sigma
  } else {
    if (is.null(mu)) {
      cli::cli_abort(
        "mu must be supplied for the normal prior",
        .internal = TRUE
      )
    }
    if (is.null(sigma)) {
      if (is.null(two_sigma)) {
        cli::cli_abort(
          "Either sigma or two_sigma must be supplied for the normal prior",
          .internal = TRUE
        )
      } else {
        sigma <- two_sigma / 2
      }
    }
  }
  return(set_prior(param = param, dist = "truncated-normal", mu = mu, sigma = sigma))
}

#' Find mu and sigma of a truncated normal distribution for given quantiles
#'
#' @description This function estimates the parameters `mu` and `sigma` of a
#'   normal distribution truncated at zero, such that the 5% and 95% quantiles
#'   of the distribution match the specified target quantiles.
#'
#' @param q5_target The target quantile at 5% (lower bound).
#' @param q95_target The target quantile at 95% (upper bound).
#'
#' @return A list containing the estimated `mu` and `sigma` for the truncated
#'   normal distribution that matches the specified quantiles.
#' @keywords internal
find_mu_sigma_truncnorm <- function(q5_target, q95_target) {
  # Check if the target quantiles are valid
  if (q5_target > q95_target) {
    cli::cli_abort(paste0(
      "Provided 5% and 95% quantiles are invalid: ",
      q5_target, " - ", q95_target, "."
    ))
  }

  if (q5_target == q95_target) {
    return(list(mu = q5_target, sigma = 0))  # fixed value
  }

  # Initial guesses
  mu_start <- (q5_target + q95_target) / 2
  sigma_start <- (q95_target - q5_target) / (qnorm(0.95) - qnorm(0.05))

  # Objective function: squared error between target and actual quantiles
  objective_fn <- function(params, q5_target, q95_target) {
    mu <- params[1]
    sigma <- params[2]
    if (sigma <= 0) return(1e10)  # penalize invalid sigma

    # Compute the 5% and 95% quantiles of the truncated normal
    q5_model <- extraDistr::qtnorm(0.05, a = 0, b = Inf, mean = mu, sd = sigma)
    q95_model <- extraDistr::qtnorm(0.95, a = 0, b = Inf, mean = mu, sd = sigma)

    # Return sum of squared errors
    sum((q5_model - q5_target)^2 + (q95_model - q95_target)^2)
  }

  # Optimize to minimize the quantile error
  result <- optim(
    par = c(mu_start, sigma_start),
    fn = objective_fn,
    q5_target = q5_target,
    q95_target = q95_target,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-6),  # sigma must be positive
    control = list(fnscale = 1)
  )

  mu <- result$par[1]
  sigma <- result$par[2]

  # check if target quantiles match realized quantiles (tolerance of 1e-2)
  if (abs(extraDistr::qtnorm(0.05, mu, sigma, a = 0) - q5_target) > 1e-3 ||
      abs(extraDistr::qtnorm(0.95, mu, sigma, a = 0) - q95_target) > 1e-3) {
    cli::cli_warn(paste0(
      "90% interval of truncated normal prior could not be exactly calibrated ",
      "(specified interval: ", q5_target, " - ", q95_target, ", ",
      "realized interval: ",
      round(extraDistr::qtnorm(0.05, mu, sigma, a = 0),3),
      " - ",
      round(extraDistr::qtnorm(0.95, mu, sigma, a = 0),3),
      "). ",
      "This typically happens when the specified quantiles are rather extreme."
    ))
  }

  return(list(mu = mu, sigma = sigma))
}

#' Provide initialization value for a parameter based on the supplied prior with location and scale
#'
#' @description Initialization using the prior is often better than initializing
#'   with zero (and if the parameter is strictly positive, zero is not possible
#'   at all). This function provides as init value the location of the prior
#'   plus 1/4 of the scale. This ensure a positive init even if the
#'   mean is zero (useful for truncated normal priors for example.)
#'
#' @param prior Prior for parameter as provided by [set_prior()]. Should be a
#'   location and scale prior (first element location, second element scale).
#'
#' @details If the provided prior has zero variance, it is assumed that the
#'   parameter will not be sampled and an empty init is returned.
#'
#' @return Location of the prior plus 1/4 of the scale.
#' @keywords internal
init_from_location_scale_prior <- function(prior, enforce_positive = FALSE) {
  prior_data_select <- !stringr::str_detect(names(prior), pattern = "_text$")
  if (length(which(prior_data_select)) != 1) {
    cli::cli_warn(paste(
      "Could not init from location scale prior:",
      "Non-ambiguous prior format. Using 1e-2 as fallback."
    ), .internal = TRUE)
  }
  prior_data <- prior[prior_data_select][[1]]
  if (prior_data[2] > 0) {
    init <- prior_data[1] + prior_data[2]/4
    if (enforce_positive && init < 0) {
      init <- prior_data[2]/4
    }
    return(init)
  } else {
    return(numeric(0))
  }
}

add_dummy_data <- function(modeldata, dummies) {
  for (dat in dummies) {
    modeldata[[dat]] <- numeric(0)
  }
  return(modeldata)
}

add_dummy_inits <- function(modeldata, dummies) {
  for (param in dummies) {
    modeldata$.init[[param]] <- numeric(0)
  }
  return(modeldata)
}

#' Print a `modeldata` object
#'
#' @param type If "structure", the model structure is printed. If "data", the
#'   modeldata content (data, inits, metainfo) is printed.
#' @export
#' @keywords internal
print.modeldata <- function(x, type = "structure", ...) {
  if (type == "structure") {
    print(x$.str)
  } else if (type == "data") {
    print.default(x[!names(x) %in% c(
      ".str", ".checks", ".spec"
      )])
  } else {
    print.default(x)
  }
}

#' Print a `modelstructure` object
#'
#' @export
#' @keywords internal
print.modelstructure <- function(x, ...) {
  output <- sapply(names(x), function(module) {
    if (length(x[[module]]) > 0) {
      comp_output <- sapply(names(x[[module]]), function(component) {
        if (length(x[[module]][[component]][[1]]) > 0) {
          a <- x[[module]][[component]][[1]]
          details <- paste0(
            " (", paste(paste0(names(a), " = ", a), collapse = ", "), ")"
            )
          paste0(" |- ", names(x[[module]][[component]])[1], details)
        } else {
          paste0(" |- ", names(x[[module]][[component]])[1])
        }
      })
      paste0(module, "\n", paste(comp_output, collapse = "\n"))
    } else {
      return(NULL)
    }
  })
  cat(paste(output[!sapply(output,is.null)], collapse = "\n\n"))
}
