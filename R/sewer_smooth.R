#' @keywords internal
configure_R_model <- function(name_approach, model_id, use_ets, use_bs,
                              use_bs2, use_scp, modeldata) {
  modeldata$.metainfo$R_estimate_approach <- name_approach
  modeldata$R_model <- model_id
  modeldata$R_use_ets <- use_ets
  modeldata$R_use_bs <- use_bs
  modeldata$R_use_bs2 <- use_bs2
  modeldata$R_use_scp <- use_scp
  return(modeldata)
}

#' @keywords internal
add_sparse_matrix <- function(M, name, modeldata){
  modeldata[[paste0(name, "_ncol")]] <- ncol(M)
  M_sparse <- suppressMessages(rstan::extract_sparse_parts(M))
  modeldata[[paste0(name, "_n_w")]] <- length(M_sparse$w)
  modeldata[[paste0(name, "_w")]] <- M_sparse$w
  modeldata[[paste0(name, "_v")]] <- M_sparse$v
  modeldata[[paste0(name, "_u")]] <- M_sparse$u
  return(modeldata)
}

#' @keywords internal
use_basis_splines <- function(spline_length, knots, degree, modeldata) {
  modeldata$bs_length <- spline_length
  modeldata$bs_dist <- median(diff(knots$interior))
  B <- splines::bs(
      1:spline_length,
      knots = knots$interior,
      degree = degree,
      intercept = FALSE,
      Boundary.knots = knots$boundary
    )
  modeldata$.metainfo$bs_matrix <- B
  modeldata <- add_sparse_matrix(B, "bs", modeldata)
  modeldata$.init$bs_coeff_noise_raw <- rep(0, modeldata$bs_ncol - 1)
  return(modeldata)
}

add_dummies_basis_splines <- function(modeldata) {
  modeldata <- add_dummy_data(modeldata, c(
    "bs_ncol", "bs_n_w", "bs_w", "bs_v", "bs_u", "bs_coeff_noise_raw"
  ))
  modeldata <- add_dummy_inits(modeldata, c(
    "bs_coeff_noise_raw"
  ))
  modeldata$bs_length <- 0
  modeldata$bs_dist <- 0
  return(modeldata)
}

#' @keywords internal
use_basis_splines2 <- function(spline_length, knots, degree, modeldata) {
  modeldata$bs2_length <- spline_length
  modeldata$bs2_dist <- median(diff(knots$interior))
  B <- splines::bs(
    1:spline_length,
    knots = knots$interior,
    degree = degree,
    intercept = FALSE,
    Boundary.knots = knots$boundary
  )
  modeldata$.metainfo$bs2_matrix <- B
  modeldata <- add_sparse_matrix(B, "bs2", modeldata)

  modeldata$.init$bs2_coeff_noise_raw <- rep(0, modeldata$bs2_ncol - 1)
  return(modeldata)
}

add_dummies_basis_splines2 <- function(modeldata) {
  modeldata <- add_dummy_data(modeldata, c(
    "bs2_ncol", "bs2_n_w", "bs2_w", "bs2_v", "bs2_u", "bs2_coeff_noise_raw"
  ))
  modeldata <- add_dummy_inits(modeldata, c(
    "bs2_coeff_noise_raw"
  ))
  modeldata$bs2_length <- 0
  modeldata$bs2_dist <- 0
  return(modeldata)
}

#' @title Solve Dirichlet Parameters for Soft Changepoint Model
#'
#' @description This function solves for the parameters p_i of a Dirichlet
#'   distribution given a total sum alpha and number of parameters n.
#'
#' @param alpha Total sum of all p_i
#' @param n Number of parameters to solve for
#' @param min_break_dist Minimum distance between changepoints. For first
#'   segment, can be left empty/set to 1.
#' @param prior_previous Vector with alpha values of prior for previous segment.
#'   Required if prior_previous>1.
#'
#' @return A vector of p_i values
#' @keywords internal
changepoint_dirichlet_prior <- function(alpha, n, min_break_dist = 1, prior_previous = NULL) {
  p <- numeric(n)
  cumulative_sum <- 0

  # Stepwise solve
  for (i in 1:(n-1)) {
    target <- i / n

    # Define function to find root for p_i
    f <- function(pi) {
      shape1 <- cumulative_sum + pi
      shape2 <- alpha - cumulative_sum - pi
      Fnow = pbeta(0.5, shape2, shape1)
      if (min_break_dist > 1 && i < min_break_dist) {
        if (is.null(prior_previous)) {
          cli::cli_abort(paste(
            "You provided a min_break_dist>1",
            "but no prior for the previous segment."
          ))
        }
        Fprevious <- pbeta(0.5, sum(prior_previous) - sum(prior_previous[1:(n - min_break_dist + i)]), sum(prior_previous[1:(n - min_break_dist + i)]))
      } else {
        Fprevious <- 1
      }
      return(Fnow * Fprevious - target)
    }

    # Use uniroot to solve for p_i in valid range
    lower <- 1e-6
    upper <- alpha - cumulative_sum - 1e-6
    result <- uniroot(f, lower = lower, upper = upper)
    p[i] <- result$root
    cumulative_sum <- cumulative_sum + p[i]
  }
  # Set last p_n to use remaining sum
  p[n] <- alpha - cumulative_sum
  return(p)
}


#' @keywords internal
use_soft_changepoints <- function(scp_length, last_knot, distance, min_distance,
                                 min_distance_tolerance, strictness_tol_k,
                                 sharpness_boltzmann, strictness_alpha,
                                 modeldata) {
  knots <- rev(seq(last_knot - distance, 1, by = -distance))
  modeldata$.metainfo$scp_knots <- knots
  modeldata$scp_break_dist <- distance
  modeldata$scp_min_dist <- min_distance
  modeldata$scp_skip_tolerance <- min_distance_tolerance
  modeldata$scp_skip_tolerance_k <- strictness_tol_k
  modeldata$scp_boltzmann_sharpness <- sharpness_boltzmann

  modeldata$scp_alpha_base <- changepoint_dirichlet_prior(
    alpha = strictness_alpha * distance,
    n = distance
    )
  modeldata$scp_alpha_adjusted <- changepoint_dirichlet_prior(
    alpha = strictness_alpha * distance,
    n = distance,
    min_break_dist = min_distance,
    prior_previous = modeldata$scp_alpha_base
  )

  modeldata$scp_n_knots <- length(knots)
  modeldata$scp_length_intercept <- knots[1] - 1
  modeldata$scp_length <- scp_length

  modeldata$.init$scp_noise <- rep(0, length(knots))
  modeldata$.init$scp_break_delays_raw <- lapply(1:length(knots),function(x) {
    rep(1e-4, distance-1)
    })
  modeldata$.init$scp_sd <- rep(1e-4, length(knots))
  return(modeldata)
}

add_dummies_soft_changepoints <- function(modeldata) {
  modeldata <- add_dummy_data(modeldata, c(
    "scp_n_knots", "scp_break_dist", "scp_min_dist",
    "scp_length_intercept", "scp_k", "scp_alpha",
    "scp_skip_tolerance", "scp_skip_tolerance_k",
    "scp_boltzmann_sharpness", "scp_alpha_base", "scp_alpha_adjusted"
  ))
  modeldata$scp_length <- 0
  modeldata <- add_dummy_inits(modeldata, c(
    "scp_noise", "scp_sd"
  ))
  modeldata$.init$scp_break_delays <- matrix(1)
  modeldata$.init$scp_break_delays_raw <- matrix(1)
  return(modeldata)
}

add_dummies_exponential_smoothing <- function(modeldata) {
  modeldata <- add_dummy_data(modeldata, c(
    "ets_trend_start_prior", "ets_diff", "ets_noncentered",
    "ets_alpha_prior", "ets_beta_prior", "ets_phi_prior"
  ))
  modeldata <- add_dummy_inits(modeldata, c(
    "ets_trend_start", "ets_noise", "ets_alpha", "ets_beta", "ets_phi"
  ))
  modeldata$ets_length <- 0
  return(modeldata)
}

add_dummies_smooth_derivative <- function(modeldata) {
  modeldata <- add_dummy_inits(modeldata, c(
    "bs_coeff_noise_sharp"
  ))
  modeldata$R_use_bs_sharp <- FALSE
  return(modeldata)
}

add_dummies_R_vari_selection <- function(modeldata) {
  R_vari_sel_dummies <- c(
    "R_vari_sel_ncol", "R_vari_sel_n_w", "R_vari_sel_w",
    "R_vari_sel_v", "R_vari_sel_u"
  )
  R_vari_sel_local_dummies <- c(
    "R_vari_sel_local_ncol", "R_vari_sel_local_n_w", "R_vari_sel_local_w",
    "R_vari_sel_local_v", "R_vari_sel_local_u"
  )
  modeldata <- add_dummy_data(modeldata, c(
    R_vari_sel_dummies, R_vari_sel_local_dummies
  ))
  return(modeldata)
}

add_dummies_R_vari <- function(modeldata, keep_baseline = FALSE) {
  R_vari_dummies <- c(
    "R_vari_ncol", "R_vari_n_w", "R_vari_w", "R_vari_v", "R_vari_u"
  )
  R_vari_dummies_inits <- c("R_sd_changepoints")
  if (!keep_baseline) {
    R_vari_dummies <- c(R_vari_dummies, "R_sd_baseline_prior")
    R_vari_dummies_inits <- c(R_vari_dummies_inits, "R_sd_baseline")
  }
  modeldata <- add_dummy_data(modeldata, c(R_vari_dummies))
  modeldata <- add_dummy_inits(modeldata, c(R_vari_dummies_inits))
  return(modeldata)
}

