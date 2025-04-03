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
  return(modeldata)
}

#' @keywords internal
use_basis_splines2 <- function(spline_length, knots, degree, modeldata) {
  modeldata$bs2_length <- spline_length
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
  return(modeldata)
}

#' @keywords internal
use_soft_changepoints <- function(scp_length, knots, distance, min_distance,
                                 min_distance_tolerance, strictness_tol_k,
                                 strictness_k, strictness_alpha,
                                 modeldata) {
  modeldata$scp_break_dist <- distance
  modeldata$scp_min_dist <- min_distance
  modeldata$scp_skip_tolerance <- min_distance_tolerance
  modeldata$scp_skip_tolerance_k <- strictness_tol_k
  modeldata$scp_k <- strictness_k
  modeldata$scp_alpha <- strictness_alpha

  modeldata$scp_n_knots <- length(knots)
  modeldata$scp_length_intercept <- knots[1] - 1
  modeldata$scp_length <- scp_length

  modeldata$.init$scp_noise <- rep(0, length(knots))
  modeldata$.init$scp_break_delays <- lapply(1:length(knots),function(x) {
    rep(1/distance,distance)
    })
  modeldata$.init$scp_sd <- rep(1e-4, length(knots))
  return(modeldata)
}

add_dummies_soft_changepoints <- function(modeldata) {
  modeldata <- add_dummy_data(modeldata, c(
    "scp_n_knots", "scp_break_dist", "scp_min_dist",
    "scp_length_intercept", "scp_k", "scp_alpha",
    "scp_skip_tolerance", "scp_skip_tolerance_k"
  ))
  modeldata$scp_length <- 0
  modeldata <- add_dummy_inits(modeldata, c(
    "scp_noise", "scp_sd"
  ))
  modeldata$.init$scp_break_delays <- matrix(1)
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

add_dummies_R_vari <- function(modeldata) {
  modeldata <- add_dummy_data(modeldata, c(
    "R_sd_baseline_prior",
    "R_vari_ncol", "R_vari_n_w", "R_vari_w", "R_vari_v", "R_vari_u"
  ))
  modeldata <- add_dummy_inits(modeldata, c(
    "R_sd_baseline", "R_sd_changepoints"
  ))
  return(modeldata)
}

