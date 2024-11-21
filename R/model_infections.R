#' Model the infection process
#'
#' @description This module function is used to specify the components of the
#'  `infections` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param generation_dist Generation time distribution. The intrinsic
#'   distribution of the time between infection of a primary case and infection
#'   of its secondary cases. Modeling options:
#' `r component_functions_("generation_dist")`
#' @param R Effective reproduction number over time. This is the main parameter
#'   of interest estimated by `EpiSewer`. `R` is smoothed using a time series
#'   smoothing prior. Currently supported are: random walk (rw), exponential
#'   smoothing (ets), and smoothing splines. Modeling options:
#' `r component_functions_("R")`
#' @param seeding Seeding of initial infections. The renewal model used by
#'   `EpiSewer` requires a seeding phase of the length of the maximum generation
#'   time. For these initial infections, a simple seeding model instead of the
#'   renewal model must be used. Modeling options:
#' `r component_functions_("seeding")`
#' @param infection_noise Noise in the infection process. `EpiSewer` implements
#' a stochastic infection model, i.e. allows for variation in the number of new
#' infections generated at each time step. This accounts for stochastic
#' uncertainty in the infection process and often speeds up model fitting.
#' Modeling options:
#' `r component_functions_("infection_noise")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `infections` module.
#' @export
#' @family {module functions}
model_infections <- function(
    generation_dist = generation_dist_assume(),
    R = R_estimate_splines(),
    seeding = seeding_estimate_rw(),
    infection_noise = infection_noise_estimate()) {
  verify_is_modeldata(generation_dist, "generation_dist")
  verify_is_modeldata(R, "R")
  verify_is_modeldata(seeding, "seeding")
  verify_is_modeldata(infection_noise, "infection_noise")
  return(modeldata_combine(generation_dist, R, seeding, infection_noise))
}

#' Assume a generation time distribution
#'
#' @description This option assumes a fixed generation time distribution for the
#'   renewal model in `EpiSewer`.
#'
#' @details The generation time distribution here refers to the intrinsic
#'   distribution of the time between infection of a primary case and infection
#'   of its secondary cases. It is disease-specific and typically obtained from
#'   literature.
#'
#' @param generation_dist A numeric vector representing a discrete generation
#'   time distribution, starting with the probability for a generation time of 1
#'   day, 2 days, 3 days, and so on (a generation time of 0 days is excluded).
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @seealso Helpers to discretize continuous probability distributions:
#'   [get_discrete_gamma()],
#'   [get_discrete_gamma_shifted()],
#'   [get_discrete_lognormal()]
generation_dist_assume <-
  function(generation_dist = NULL, modeldata = modeldata_init()) {
    modeldata <- tbp("generation_dist_assume",
      {
        modeldata$G <- length(generation_dist)
        generation_dist <- check_dist(
          generation_dist,
          "generation time distribution"
        )
        modeldata$generation_dist <- generation_dist
      },
      required_assumptions = "generation_dist",
      modeldata = modeldata
    )

    modeldata$.str$infections[["generation_dist"]] <- list(
      generation_dist_assume = c()
    )

    return(modeldata)
  }

#'Estimate Rt via exponential smoothing
#'
#'@description This option estimates the effective reproduction number over time
#'  using exponential smoothing. It implements Holt's linear trend method with
#'  dampening through an innovations state space model with a level, trend, and
#'  dampening component.
#'
#'@param level_prior_mu Prior (mean) on the initial level of Rt.
#'@param level_prior_sigma Prior (standard deviation) on the initial level of
#'  Rt.
#'@param trend_prior_mu Prior (mean) on the initial trend of Rt.
#'@param trend_prior_sigma Prior (standard deviation) on the initial trend of
#'  Rt.
#'@param sd_base_prior_sd Prior (standard deviation) on the baseline standard
#'  deviation of the innovations. We here use a half-normal prior, i.e.
#'  `sd_base_prior_sd` is the only parameter to be specified for this prior.
#'  Please note that for consistency, the overall standard deviation of
#'  innovations will always be the baseline plus an additive component from
#'  `sd_change_prior` even if no changepoints are modeled (see below).
#'@param sd_change_distance Distance between changepoints used to model
#'  additional variation in Rt. The default change point distance is 4 weeks.
#'  Very short changepoint distances must be chosen with care, as they can make
#'  the Rt time series too flexible. If set to zero, no change points are
#'  modeled.
#'@param sd_change_prior_shape Exponential-Gamma prior (shape) on standard
#'  deviation additional to baseline. This prior describes the distribution of
#'  the standard deviation of Rt over time. EpiSewer will estimate a baseline
#'  standard deviation (see `sd_base_prior_sd`), and model additional variation
#'  on top of the baseline using a changepoint model. Please see the details for
#'  more explanation.
#'@param sd_change_prior_rate Exponential-Gamma prior (rate) on standard
#'  deviation additional to baseline. See `sd_change_prior_shape` and the
#'  details for more explanation.
#'@param smooth_prior_mu Prior (mean) on the smoothing parameter. Must be
#'  between 0 and 1.
#'@param smooth_prior_sigma Prior (standard deviation) on the smoothing
#'  parameter. If this is set to zero, the smoothing parameter will be fixed to
#'  `smooth_prior_mu` and not estimated. If positive, a beta prior with the
#'  corresponding mean and standard deviation is used.
#'@param trend_smooth_prior_mu Prior (mean) on the trend smoothing parameter.
#'  Must be between 0 and 1.
#'@param trend_smooth_prior_sigma Prior (standard deviation) on the trend
#'  smoothing parameter. If this is set to zero, the trend smoothing parameter
#'  will be fixed to `trend_smooth_prior_mu` and not estimated. If positive, a
#'  beta prior with the corresponding mean and standard deviation is used.
#'@param dampen_prior_mu Prior (mean) on the dampening parameter. Must be
#'  between 0 and 1.
#'@param dampen_prior_sigma Prior (standard deviation) on the dampening
#'  parameter. If this is set to zero, the dampening parameter will be fixed to
#'  `dampen_prior_mu` and not estimated. If positive, a beta prior with the
#'  corresponding mean and standard deviation is used.
#'@param differenced If `FALSE` (default), exponential smoothing is applied to
#'  the absolute Rt time series. If `TRUE`, it is instead applied to the
#'  differenced time series. This makes the level become the trend, and the
#'  trend become the curvature.
#'@param noncentered If `TRUE` (default), a non-centered parameterization is
#'  used to model the innovations in the state space process (for better
#'  sampling efficiency).
#'
#'@inheritParams R_estimate_splines
#'
#'@details The innovations state space model consists of three components: a
#'  level, a trend, and a dampening component.
#' - The level is smoothed based on the levels from earlier time steps,
#'  with exponentially decaying weights, as controlled by a smoothing parameter
#'  (often called alpha). Note that *smaller* values of `alpha` indicate
#'   *stronger* smoothing. In particular, `alpha = 1` means that only the last
#'  level is used.
#' - The trend is smoothed based on the trends from earlier time steps,
#'  with exponentially decaying weights, as controlled by a trend smoothing
#'  parameter (often called beta). Note that *smaller* values of `beta` indicate
#'  *stronger* smoothing. In particular, `beta = 1` means that only the last
#'  trend is used.
#' - The dampening determines how long a previous trend continues into the
#'  future before it levels of to a stationary time series. The strength of
#'  dampening is controlled by a dampening parameter (often called phi). Note
#'  that *smaller* values of `phi` indicate *stronger* dampening. In particular,
#'  `phi = 1` means no dampening. Values below `phi = 0.8` are seldom in
#'  practice as the dampening becomes very strong.
#'
#'@details Often, `alpha`, `beta`, and `phi` are jointly unidentifiable. It may
#'  therefore be necessary to fix at least one of the parameters (typically
#'  `phi`) or supply strong priors.
#'
#'@details Note that the smoothness of retrospective Rt estimates is often more
#'  influenced by the prior on the standard deviation of innovations than the
#'  smoothing and trend smoothing parameters. The smoothing parameters mostly
#'  have an influence on estimates close to the present / date of estimation,
#'  when limited data signal is available. Here, the standard deviation of the
#'  innovations influences how uncertain Rt estimates are close to the present.
#'
#'@details The variability of Rt can change over time. For example, during the
#'  height of an epidemic wave, countermeasures may lead to much faster changes
#'  in Rt than observable at other times (baseline). This potential additional
#'  variability is accounted for using change points placed at regular
#'  intervals. The additional standard deviation of the state space model
#'  innovations on top of the baseline then evolves linearly between the change
#'  points. The additional variation defined at the changepoints is modeled as
#'  independently distributed and following a Lomax distribution, also known as
#'  Exponential-Gamma (EG) distribution. This is an exponential distribution
#'  where the rate is Gamma distributed. The prior `sd_change_prior` defines the
#'  shape and rate of this Gamma distribution. The distribution has a
#'  strong peak towards zero and a long tail. This regularizes the estimated
#'  deviations from the baseline standard deviation - most deviations are small,
#'  but during special time periods, the deviation might also be larger.
#'
#'@details The priors of this component have the following functional form:
#' - initial level of Rt: `Normal`
#' - initial trend of Rt: `Normal`
#' - baseline standard deviation of innovations: `Half-normal`
#' - additional standard deviation at changepoints: `Exponential-Gamma`
#' - smoothing parameter: `Beta`
#' - trend smoothing parameter: `Beta`
#' - dampening parameter: `Beta`
#'
#'@inheritParams template_model_helpers
#'@inherit modeldata_init return
#'@export
#'@family {Rt models}
R_estimate_ets <- function(
    level_prior_mu = 1,
    level_prior_sigma = 0.8,
    trend_prior_mu = 0,
    trend_prior_sigma = 0.1,
    sd_base_prior_sd = 0.025,
    sd_change_prior_shape = 0.5,
    sd_change_prior_rate = 1e-4,
    sd_change_distance = 7*26,
    link = "inv_softplus",
    R_max = 6,
    smooth_prior_mu = 0.5,
    smooth_prior_sigma = 0.05,
    trend_smooth_prior_mu = 0.5,
    trend_smooth_prior_sigma = 0.05,
    dampen_prior_mu = 0.9,
    dampen_prior_sigma = 0,
    differenced = FALSE,
    noncentered = TRUE,
    modeldata = modeldata_init()) {
  modeldata$.metainfo$R_estimate_approach <- "ets"
  modeldata$R_model <- 0

  modeldata$R_level_start_prior <- set_prior(
    "R_level_start", "normal",
    mu = level_prior_mu, sigma = level_prior_sigma
  )

  modeldata$R_trend_start_prior <- set_prior(
    "R_trend_start", "normal",
    mu = trend_prior_mu, sigma = trend_prior_sigma
  )

  modeldata$R_sd_baseline_prior <- set_prior("R_sd_baseline",
    "half-normal",
    sigma = sd_base_prior_sd
  )
  modeldata$R_sd_change_prior <- set_prior("R_sd_change",
    "lomax",
    shape = sd_change_prior_shape,
    rate = sd_change_prior_rate
  )

  modeldata$.init$R_level_start <-
    modeldata$R_level_start_prior$R_level_start_prior[1]
  modeldata$.init$R_trend_start <- 1e-4

  modeldata <- tbc(
    "R_ets_noise",
    {
      modeldata$.init$R_noise <- rep(0, modeldata$.metainfo$length_R - 1)
      modeldata <- add_R_variability(
        length_R = modeldata$.metainfo$length_R,
        length_seeding = modeldata$.metainfo$length_seeding,
        length_partial = modeldata$.metainfo$partial_window,
        length_partial_gen = modeldata$.metainfo$partial_generation,
        h = modeldata$.metainfo$forecast_horizon,
        changepoint_dist = sd_changepoint_dist,
        modeldata = modeldata
      )
      modeldata$.init$R_sd_baseline <- 1e-2
      if (modeldata$R_vari_n_basis > 2) {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_n_basis - 2)
      } else {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_n_basis)
      }
    },
    required = c(
      ".metainfo$length_R",
      ".metainfo$length_seeding",
      ".metainfo$partial_window",
      ".metainfo$partial_generation",
      ".metainfo$forecast_horizon"
    ),
    modeldata = modeldata
  )

  check_beta_alternative(smooth_prior_mu, smooth_prior_sigma)
  modeldata$ets_alpha_prior <- set_prior("ets_alpha", "beta",
    mu = smooth_prior_mu, sigma = smooth_prior_sigma
  )
  modeldata$.init$ets_alpha <- init_from_location_scale_prior(
    modeldata$ets_alpha_prior
  )

  check_beta_alternative(trend_smooth_prior_mu, trend_smooth_prior_sigma)
  modeldata$ets_beta_prior <- set_prior("ets_beta", "beta",
    mu = trend_smooth_prior_mu, sigma = trend_smooth_prior_sigma
  )
  modeldata$.init$ets_beta <- init_from_location_scale_prior(
    modeldata$ets_beta_prior
  )

  check_beta_alternative(dampen_prior_mu, dampen_prior_sigma)
  modeldata$ets_phi_prior <- set_prior("ets_phi", "beta",
    mu = dampen_prior_mu, beta = dampen_prior_sigma
  )
  modeldata$.init$ets_phi <- init_from_location_scale_prior(
    modeldata$ets_phi_prior
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  modeldata$ets_diff <- differenced
  modeldata$ets_noncentered <- noncentered

  modeldata <- add_dummy_data(modeldata, c(
    "bs_n_basis", "bs_n_w", "bs_w", "bs_v", "bs_u",
    "bs_coeff_ar_start_prior",
    "R_vari_sel_ncol", "R_vari_sel_n_w", "R_vari_sel_w", "R_vari_sel_v",
    "R_vari_sel_u"
  ))

  modeldata <- add_dummy_inits(modeldata, c(
    "bs_coeff_ar_start", "bs_coeff_ar_sd", "bs_coeff_noise_raw"
  ))

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_ets = c()
  )

  return(modeldata)
}

#'Estimate Rt via a random walk
#'
#'@description This option estimates the effective reproduction number over time
#'  using a random walk.
#'
#'@param intercept_prior_mu Prior (mean) on the intercept of the random walk.
#'@param intercept_prior_sigma Prior (standard deviation) on the intercept of
#'  the random walk.
#'@param sd_base_prior_sd Prior (standard deviation) on the baseline standard
#'  deviation of the innovations. We here use a half-normal prior, i.e.
#'  `sd_base_prior_sd` is the only parameter to be specified for this prior.
#'  Please note that for consistency, the overall standard deviation of
#'  innovations will always be the baseline plus an additive component from
#'  `sd_change_prior` even if no changepoints are modeled (see below).
#'@param sd_change_distance Distance between changepoints used to model
#'  additional variation in Rt. The default change point distance is 4 weeks.
#'  Very short changepoint distances must be chosen with care, as they can make
#'  the Rt time series too flexible. If set to zero, no change points are
#'  modeled.
#'@param sd_change_prior_shape Exponential-Gamma prior (shape) on standard
#'  deviation additional to baseline. This prior describes the distribution of
#'  the standard deviation of Rt over time. EpiSewer will estimate a baseline
#'  standard deviation (see `sd_base_prior_sd`), and model additional variation
#'  on top of the baseline using a changepoint model. Please see the details for
#'  more explanation.
#'@param sd_change_prior_rate Exponential-Gamma prior (rate) on standard
#'  deviation additional to baseline. See `sd_change_prior_shape` and the
#'  details for more explanation.
#'@param differenced If `FALSE` (default), the random walk is applied to the
#'  absolute Rt time series. If `TRUE`, it is instead applied to the differenced
#'  time series, i.e. now the trend is modeled as a random walk.
#'@param noncentered If `TRUE` (default), a non-centered parameterization is
#'  used to model the innovations of the random walk (for better sampling
#'  efficiency).
#'@inheritParams R_estimate_ets
#'
#'@details The smoothness of Rt estimates is influenced by the prior on the
#'  standard deviation of the random walk. It also influences the uncertainty of
#'  Rt estimates towards the present / date of estimation, when limited data
#'  signal is available. The prior on the intercept of the random walk should
#'  reflect your expectation of Rt at the beginning of the time series. If
#'  estimating from the start of an epidemic, you might want to use a prior with
#'  mean > 1 for the intercept.
#'
#'@details The variability of Rt can change over time. For example, during the
#'  height of an epidemic wave, countermeasures may lead to much faster changes
#'  in Rt than observable at other times (baseline). This potential additional
#'  variability is accounted for using change points placed at regular
#'  intervals. The additional standard deviation of the state space model
#'  innovations on top of the baseline then evolves linearly between the change
#'  points. The additional variation defined at the changepoints is modeled as
#'  independently distributed and following a Lomax distribution, also known as
#'  Exponential-Gamma (EG) distribution. This is an exponential distribution
#'  where the rate is Gamma distributed. The prior `sd_change_prior` defines the
#'  shape and rate of this Gamma distribution. The distribution has a
#'  strong peak towards zero and a long tail. This regularizes the estimated
#'  deviations from the baseline standard deviation - most deviations are small,
#'  but during special time periods, the deviation might also be larger.
#'
#'@details The priors of this component have the following functional form:
#' - intercept of the random walk: `Normal`
#' - baseline standard deviation of the random walk: `Half-normal`
#' - additional standard deviation at changepoints: `Exponential-Gamma`
#'
#'@inheritParams template_model_helpers
#'@inherit modeldata_init return
#'@export
#'@family {Rt models}
R_estimate_rw <- function(
    intercept_prior_mu = 1,
    intercept_prior_sigma = 0.8,
    sd_base_prior_sd = 0.025,
    sd_change_prior_shape = 0.5,
    sd_change_prior_rate = 1e-4,
    sd_change_distance = 7*26,
    link = "inv_softplus",
    R_max = 6,
    differenced = FALSE,
    noncentered = TRUE,
    modeldata = modeldata_init()) {
  modeldata <- R_estimate_ets(
    level_prior_mu = intercept_prior_mu,
    level_prior_sigma = intercept_prior_sigma,
    sd_base_prior_sd = sd_base_prior_sd,
    sd_change_prior_shape = sd_change_prior_shape,
    sd_change_prior_rate = sd_change_prior_rate,
    sd_change_distance = sd_change_distance,
    link = link,
    R_max = R_max,
    smooth_prior_mu = 1,
    smooth_prior_sigma = 0,
    trend_smooth_prior_mu = 0,
    trend_smooth_prior_sigma = 0,
    dampen_prior_mu = 0,
    dampen_prior_sigma = 0,
    differenced = differenced,
    noncentered = noncentered,
    modeldata = modeldata
  )

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_rw = c()
  )

  return(modeldata)
}

#'Estimate Rt via smoothing splines
#'
#'@description This option estimates the effective reproduction number using
#'  penalized B-splines.
#'
#'@param knot_distance Distance between spline breakpoints (knots) in days
#'  (default is 7, i.e. one knot each week). Placing knots further apart
#'  increases the smoothness of Rt estimates and can speed up model fitting. The
#'  Rt time series remains surprisingly flexible even at larger knot distances,
#'  but placing knots too far apart can lead to inaccurate estimates.
#'@param spline_degree Degree of the spline polynomials (default is 3 for cubic
#'  splines).
#'@param R_intercept_prior_mu Prior (mean) on the intercept of the random
#'  walk over spline coefficients.
#'@param R_intercept_prior_sigma Prior (standard deviation) on the intercept
#'  of the random walk over spline coefficients.
#'@param R_sd_local_prior_sd Prior (standard deviation) on the baseline standard
#'  deviation of the random walk over spline coefficients. We here use a
#'  half-normal prior, i.e. `sd_base_prior_sd` is the only parameter to be
#'  specified for this prior. Please note that for consistency, the overall
#'  standard deviation of the random walk over spline coefficients will always
#'  be the baseline plus an additive component from `sd_change_prior` even if no
#'  changepoints are modeled (see below).
#'@param coef_sd_change_distance Distance between changepoints used to model
#'  additional variation in Rt. The default change point distance is 4 weeks.
#'  Very short changepoint distances must be chosen with care, as they can make
#'  the Rt time series too flexible. If set to zero, no change points are
#'  modeled.
#'@param coef_sd_change_prior_shape Exponential-Gamma prior (shape) on standard
#'  deviation additional to baseline. This prior describes the distribution of
#'  the standard deviation of the random walk over spline coefficients. EpiSewer
#'  will estimate a baseline standard deviation (see `sd_base_prior_sd`), and
#'  model additional variation on top of the baseline using a changepoint model.
#'  Please see the details for more explanation.
#'@param coef_sd_change_prior_rate Exponential-Gamma prior (rate) on standard
#'  deviation additional to baseline. See `sd_change_prior_shape` and the
#'  details for more explanation.
#'@param link Link function. Currently supported are `inv_softplus` (default)
#'  and `scaled_logit`. Both of these links are configured to behave
#'  approximately like the identity function around R=1, but become increasingly
#'  non-linear below (and in the case of `scaled_logit` also above) R=1.
#'@param R_max If `link=scaled_logit` is used, a maximum reproduction number
#'  must be assumed. This should be higher than any realistic R value for the
#'  modeled pathogen. Default is 6.
#'
#'@details `EpiSewer` uses a random walk on the B-spline coefficients for
#'  regularization. This allows to use small knot distances without obtaining
#'  extremely wiggly Rt estimates.
#'   - The prior on the random walk intercept should reflect
#'  your expectation of Rt at the beginning of the time series. If estimating
#'  from the start of an epidemic, you might want to use a prior with mean
#'  larger than 1 for the intercept.
#'   - The prior on the baseline standard deviation should be interpreted in terms
#'  of daily additive changes (this is accurate around Rt=1, and becomes less
#'  accurate as Rt approaches 0 or its upper bound as defined by the `link`
#'  function). For example, a baseline half-normal prior with sd=0.05 allows a
#'  daily standard deviation between 0 and 0.1. A daily standard deviation of
#'  0.1 in turn roughly allows the spline coefficients to change by ±0.2 (using
#'  the 2 sigma rule) each day. The daily standard deviation is summed up over
#'  the days between two knots to get the actual standard deviation of the
#'  coefficients. This way, the prior is independent of the chosen
#'  `knot_distance`. For example, if `knot_distance` is 7 days, and a constant
#'  daily standard deviation of 0.1 is estimated, the coefficients of two
#'  adjacent splines can differ by up to `0.2*sqrt(knot_distance)`, i.e. ±0.5.
#'  Note however that this difference does not directly translate into a change
#'  of Rt by ±0.5, as Rt is always the weighted sum of several basis functions
#'  at any given point. It may therefore change much more gradually, depending
#'  on the distances between knots.
#'  - The variability of Rt can change over time. For example, during the height
#'  of an epidemic wave, countermeasures may lead to much faster changes in Rt
#'  than observable at other times (baseline). This potential additional
#'  variability is accounted for using change points placed at regular
#'  intervals. The additional standard deviation of the random walk over spline
#'  coefficients (on top of the baseline) then evolves linearly between the
#'  change points. The values at the changepoints are modeled as independently
#'  distributed and following a Lomax distribution, also known as
#'  Exponential-Gamma (EG) distribution. This is an exponential distribution
#'  where the rate is Gamma distributed. The prior `sd_change_prior` defines the
#'  shape and rate of this Gamma distribution. The distribution has a
#'  strong peak towards zero and a long tail. This regularizes the estimated
#'  deviations from the baseline standard deviation - most deviations are small,
#'  but during special time periods, the deviation might also be larger.
#'
#'@details The smoothness of the Rt estimates is influenced both by the knot
#'  distance and by the daily standard deviation of the random walk on
#'  coefficients. The latter also influences the uncertainty of Rt estimates
#'  towards the present / date of estimation, when limited data signal is
#'  available. Absent sufficient data signal, Rt estimates will tend to stay at
#'  the current level (which corresponds to assuming unchanged transmission
#'  dynamics).
#'
#'@details The priors of this component have the following functional form:
#' - intercept of the random walk over spline coefficients:
#'  `Normal`
#' - standard deviation of the random walk over spline coefficients:
#'  `Truncated normal`
#'
#'@inheritParams template_model_helpers
#'@inherit modeldata_init return
#'@export
#'@family {Rt models}
R_estimate_splines <- function(
    knot_distance_global = 4*7,
    knot_distance_local = 7,
    R_intercept_prior_mu = 1,
    R_intercept_prior_sigma = 0.8,
    R_sd_local_prior_sd = 0.05,
    coef_sd_change_prior_shape = 0.5,
    coef_sd_change_prior_rate = 1e-4,
    coef_sd_change_distance = 4*7,
    spline_degree = 3,
    link = "inv_softplus",
    R_max = 6,
    modeldata = modeldata_init()) {
  modeldata$.metainfo$R_estimate_approach <- "splines"
  modeldata$R_model <- 1

  modeldata <- tbc(
    "spline_definition",
    {
      # Global spline model for Rt
      # knots_global <- place_knots(
      #   ts_length = with(
      #     modeldata$.metainfo, length_R + forecast_horizon
      #   ),
      #   knot_distance = knot_distance_global,
      #   partial_window = max(with(
      #     modeldata$.metainfo, partial_window + forecast_horizon
      #     ), knot_distance_global + modeldata$.metainfo$forecast_horizon)
      # )
      ts_length = with(
        modeldata$.metainfo, length_R + forecast_horizon
        )
      m_meta <- modeldata$.metainfo
      partial_knots_dists <- place_knots_partial_window(m_meta$partial_window, 3)
      knots_global <- list(
        interior = rev(c(
          min(ts_length,m_meta$length_R+3),
          seq(m_meta$length_R+2, m_meta$length_R, by = -1), # forecast knots
          m_meta$length_R-partial_knots_dists, # knots during partial window
          # next knot has a distance of 90% of the generation time distribution
          seq(m_meta$length_R-partial_knots_dists[3]-m_meta$partial_generation, 1, by = -knot_distance_global)
          )),
        boundary = c(0, ts_length + 1)
      )
      B_global <-
        splines::bs(
          with(modeldata$.metainfo, 1:(length_R + forecast_horizon)),
          knots = knots_global$interior,
          degree = spline_degree,
          intercept = FALSE,
          Boundary.knots = knots_global$boundary
        )
      modeldata$.metainfo$R_knots_global <- knots_global
      modeldata$.metainfo$B_global <- B_global
      modeldata$bsg_n_basis <- ncol(B_global)
      B_sparse_global <- suppressMessages(rstan::extract_sparse_parts(B_global))
      modeldata$bsg_n_w <- length(B_sparse_global$w)
      modeldata$bsg_w <- B_sparse_global$w
      modeldata$bsg_v <- B_sparse_global$v
      modeldata$bsg_u <- B_sparse_global$u

      modeldata$.init$bsg_coeff_ar_start <- 1
      modeldata$.init$bsg_coeff_noise_raw <- rep(0, modeldata$bsg_n_basis - 1)

      # Local spline model for Rt
      # knots_local <- place_knots(
      #   ts_length = with(
      #     modeldata$.metainfo, length_R + forecast_horizon
      #   ),
      #   knot_distance = knot_distance_local,
      #   partial_window = max(with(
      #     modeldata$.metainfo, partial_window + forecast_horizon
      #   ), knot_distance_local + modeldata$.metainfo$forecast_horizon)
      # )
      ts_length = with(
        modeldata$.metainfo, length_R + forecast_horizon
      )
      knots_local <- list(
        interior = rev(c(
          min(ts_length,m_meta$length_R+3),
          seq(m_meta$length_R+2, m_meta$length_R, by = -1), # forecast knots
          m_meta$length_R-partial_knots_dists, # knots during partial window
          # next knot has a distance of 90% of the generation time distribution
          seq(m_meta$length_R-partial_knots_dists[3]-m_meta$partial_generation, 1, by = -knot_distance_local)
        )),
        boundary = c(0, ts_length + 1)
      )
      B_local <-
        splines::bs(
          with(modeldata$.metainfo, 1:(length_R + forecast_horizon)),
          knots = knots_local$interior,
          degree = spline_degree,
          intercept = FALSE,
          Boundary.knots = knots_local$boundary
        )
      modeldata$.metainfo$R_knots_local <- knots_local
      modeldata$.metainfo$B_local <- B_local
      modeldata$bsc_n_basis <- ncol(B_local)
      B_sparse_local <- suppressMessages(rstan::extract_sparse_parts(B_local))
      modeldata$bsc_n_w <- length(B_sparse_local$w)
      modeldata$bsc_w <- B_sparse_local$w
      modeldata$bsc_v <- B_sparse_local$v
      modeldata$bsc_u <- B_sparse_local$u

      modeldata$.init$bsc_coeff_noise_raw <- rep(0, modeldata$bsc_n_basis - 1)

      # Changepoint model for variability of Rt
      modeldata <- add_R_variability(
        length_R = modeldata$.metainfo$length_R,
        length_seeding = modeldata$.metainfo$length_seeding,
        length_partial = modeldata$.metainfo$partial_window,
        length_partial_gen = modeldata$.metainfo$partial_generation,
        h = modeldata$.metainfo$forecast_horizon,
        changepoint_dist = coef_sd_change_distance,
        modeldata = modeldata
        )
      modeldata$.init$R_sd_baseline <- 1e-2
      if (modeldata$R_vari_n_basis > 2) {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_n_basis - 2)
      } else {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_n_basis)
      }

      # build sparse matrix that represents which days need to be summed up
      # to compute the variance between two knots
      all_positions <- c(
        rev(knots_global$interior[1] - seq(
          0,
          by = (knots_global$interior[1] - knots_global$boundary[1]),
          length.out = spline_degree
          )), knots_global$interior[-1], knots_global$boundary[2])
      #all_positions <- c(knots_global$interior, knots_global$boundary[2])
      R_vari_selection <- t(mapply(
        get_selection_vector,
        from = all_positions[-length(all_positions)] + 1,
        to = all_positions[-1],
        n = nrow(B_global)
      ))
      modeldata$R_vari_sel_ncol <- ncol(R_vari_selection)
      R_vari_selection_sparse <- suppressMessages(
        rstan::extract_sparse_parts(R_vari_selection)
        )
      modeldata$R_vari_sel_n_w <- length(R_vari_selection_sparse$w)
      modeldata$R_vari_sel_w <- R_vari_selection_sparse$w
      modeldata$R_vari_sel_v <- R_vari_selection_sparse$v
      modeldata$R_vari_sel_u <- R_vari_selection_sparse$u

      all_positions_local <- c(
        rev(knots_local$interior[1] - seq(
          0,
          by = (knots_local$interior[1] - knots_local$boundary[1]),
          length.out = spline_degree
        )), knots_local$interior[-1], knots_local$boundary[2])
      R_vari_selection_local <- t(mapply(
        get_selection_vector,
        from = all_positions_local[-length(all_positions_local)] + 1,
        to = all_positions_local[-1],
        n = nrow(B_local)
      ))
      partial_shift <- partial_knots_dists[3] + m_meta$partial_generation
      R_vari_scaling_local <- c(
        rep(1, m_meta$length_R-partial_shift),
        ((partial_shift-1):0)/partial_shift,
        rep(0,m_meta$forecast_horizon)
        )
      # multiply each row with the scaling factor
      R_vari_selection_local <- t(t(R_vari_selection_local) * R_vari_scaling_local)
      modeldata$R_vari_sel_local_ncol <- ncol(R_vari_selection_local)
      R_vari_selection_local_sparse <- suppressMessages(
        rstan::extract_sparse_parts(R_vari_selection_local)
      )
      modeldata$R_vari_sel_local_n_w <- length(R_vari_selection_local_sparse$w)
      modeldata$R_vari_sel_local_w <- R_vari_selection_local_sparse$w
      modeldata$R_vari_sel_local_v <- R_vari_selection_local_sparse$v
      modeldata$R_vari_sel_local_u <- R_vari_selection_local_sparse$u
    },
    required = c(
      ".metainfo$length_R",
      ".metainfo$length_seeding",
      ".metainfo$partial_window",
      ".metainfo$partial_generation",
      ".metainfo$forecast_horizon"
      ),
    modeldata = modeldata
  )

  modeldata$bsg_coeff_ar_start_prior <- set_prior("bsg_coeff_ar_start",
    "normal",
    mu = R_intercept_prior_mu,
    sigma = R_intercept_prior_sigma
  )
  modeldata$R_sd_baseline_prior <- set_prior("R_sd_baseline",
    "half-normal",
    sigma = R_sd_local_prior_sd
    )
  modeldata$R_sd_change_prior <- set_prior("R_sd_change",
    "lomax",
    shape = coef_sd_change_prior_shape,
    rate = coef_sd_change_prior_rate
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  modeldata <- add_dummy_data(modeldata, c(
    "R_level_start_prior", "R_trend_start_prior",
    "ets_diff", "ets_noncentered",
    "ets_alpha_prior",
    "ets_beta_prior",
    "ets_phi_prior"
    ))

  modeldata <- add_dummy_inits(modeldata, c(
    "R_level_start", "R_trend_start", "R_noise",
    "ets_alpha", "ets_beta", "ets_phi"
    ))

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_splines = c()
  )

  return(modeldata)
}

#' Estimate Rt using an approximation of the generative renewal model
#'
#' @description This option estimates the effective reproduction number Rt using
#'   an approximation of the generative renewal model. It can be considerably
#'   faster than the fully generative options, but is often less exact. It is
#'   recommended to check the results of this method against a fully generative
#'   version like [R_estimate_splines()] before using it for real-time R
#'   estimation. See the details for an explanation of the method and its
#'   limitations.
#'
#' @param inf_sd_prior_mu Prior (mean) on the standard deviation of daily
#'   changes in infections. Infection numbers are modeled on the log scale, thus
#'   the standard deviation can be roughly interpreted in terms of percentage
#'   changes. A higher standard deviation will lead to more volatile infections
#'   and more uncertainty in R estimates. The default (`inf_sd_prior_mu = 0.05`)
#'   allows daily growth/decline rates of up +-10%. We suggest to only decrease
#'   this if you are very certain that the modeled scenario has lower maximum
#'   growth/decline rates, and to increase this if you expect a more extreme
#'   scenario. Note that the smoothness of the infection curve is also
#'   influenced by the exponential smoothing parameters (`inf_smooth`,
#'   `inf_trend_smooth`). Still, `inf_sd_prior_mu` is the primary hyperparameter
#'   to adjust if you think the estimated infection trajectory is too volatile /
#'   not flexible enough.
#' @param inf_sd_prior_sigma Prior (standard deviation) on standard deviation of
#'   daily changes of infections. This reflects uncertainty around
#'   `inf_sd_prior_mu`. The default (`inf_sd_prior_mu=0.025`) allows the
#'   standard deviation to be roughly 0.5 higher or lower than
#'   `inf_sd_prior_mu`. For the default setting with `inf_sd_prior_mu = 0.05`,
#'   this means that growth/decline rates of up to +-20% are still supported by
#'   the prior.
#' @param inf_smooth Smoothing parameter (alpha) for infections. The infection
#'   time series is smoothed based on earlier time steps, with exponentially
#'   decaying weights, as controlled by the alpha parameter. The default is 0.5,
#'   *smaller* values will lead to *stronger* smoothing, and *larger* values to
#'   *less* smoothing.
#' @param inf_trend_smooth Trend smoothing parameter (beta) for infections. The
#'   exponential trend in infections is also smoothed based on the trend over
#'   earlier time steps, with exponentially decaying weights, as controlled by
#'   the beta parameter. Default is 0.5, *smaller* values will lead to
#'   *stronger* smoothing of the trend (useful for longer generation times),
#'   and *larger* values to *less* smoothing of the trend (useful for shorter
#'   generation times).
#' @param inf_trend_dampen Trend dampening parameter (phi) for infections. The
#'   exponential trend in infections is dampened over time (exponential growth
#'   will not continue for ever). This mainly influences the trend towards the
#'   present. For the default value (0.9), the trend will only roughly half as
#'   strong after two weeks. Note that increasing phi to close to 1 can reduce
#'   sampling speed.
#' @param knot_distance The infection time series is smoothed using penalized
#'   B-splines. The distance between spline breakpoints (knots) is 7 days by
#'   default. Placing knots further apart increases the smoothness of the
#'   infection time series and can speed up model fitting. Placing knots closer
#'   to each other increases the volatility of the infections. Note however that
#'   the uncertainty of R estimates produced by [R_estimate_approx()] is also
#'   sensitive to the smoothness of the infection time series, i.e. changing the
#'   knot distance can make the R estimates appear more or less uncertain. This
#'   effect is however rather small since the regularization of the spline
#'   coefficients adjusts for the knot distance.
#' @param spline_degree Degree of the spline polynomials (default is 3 for cubic
#'   splines).
#' @param R_window Smoothing window size for R estimates. Default is 1 (i.e. no
#'   smoothing). This is provided for compatibility with the EpiEstim method by
#'   Cori et al., which assumes that R stays constant over a certain time
#'   window. However, this can artificially reduce the uncertainty of the R
#'   estimates and is therefore generally not recommended. It also means that R
#'   cannot be estimated up to the present because the smoothing window must be
#'   centered.
#'
#' @details This method uses an approximation of the renewal model to estimate
#'   the effective reproduction number. Instead of generating the infections
#'   through a renewal process, the infection time series is estimated using a
#'   non-parametric smoothing prior that mimics the characteristics of a renewal
#'   process. Rt is then computed from the infection time series by applying the
#'   classical renewal equation. This still gives Bayesian estimates of Rt and
#'   the infection time series, but is less exact than a fully generative model.
#'
#' @details To smooth the infection time series, penalized B-splines similar to
#'   the ones in [R_estimate_splines()] are used, but they are applied to the
#'   expected number of infections on the log scale, not Rt. The spline
#'   coefficients are themselves smoothed using an exponential smoothing prior
#'   with a trend component. The trend component is important as it reflects the
#'   default assumption of constant transmission dynamics that leads to
#'   exponential growth in infections. This ensures consistent results of
#'   `R_estimate_approx()` also towards the present.
#'
#' @details The method can also model infection noise as specified by
#'   [infection_noise_estimate()]. To account for the autocorrelation of
#'   infection noise, a correction that applies the renewal model specifically
#'   to the noise component is used. This correction is most accurate when R is
#'   close to 1 and less accurate when it is far from 1.
#'
#' @details The main limitation of `R_estimate_approx()` is that its priors and
#'   hyperparameters are less interpretable and thus difficult to specify. In
#'   particular the assumed generation time distribution is only used for R
#'   estimation but does not inform the smoothness of the expected infection
#'   time series. Instead, the smoothness is determined by the spline settings
#'   and exponential smoothing parameters, which are difficult to translate into
#'   a generation time assumption. The best approach is therefore to compare
#'   estimates from `R_estimate_approx()` with those from another method like
#'   [R_estimate_splines()] for a pathogen of interest and to carefully adjust
#'   the smoothing hyperparameters if needed.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
R_estimate_approx <- function(
    inf_sd_prior_mu = 0.05,
    inf_sd_prior_sigma = 0.025,
    inf_smooth = 0.75,
    inf_trend_smooth = 0.75,
    inf_trend_dampen = 0.9,
    knot_distance = 7,
    spline_degree = 3,
    R_window = 1,
    modeldata = modeldata_init()) {
  modeldata$.metainfo$R_estimate_approach <- "approx"

  modeldata <- tbc(
    "spline_definition",
    {
      knots <- place_knots(
        ts_length = with(modeldata$.metainfo, length_I - length_seeding + forecast_horizon),
        knot_distance = knot_distance,
        partial_window = with(modeldata$.metainfo, partial_window + forecast_horizon)
      )
      B <-
        splines::bs(
          1:with(modeldata$.metainfo, length_I - length_seeding + forecast_horizon),
          knots = knots$interior,
          degree = spline_degree,
          intercept = FALSE,
          Boundary.knots = knots$boundary
        )
      modeldata$.metainfo$R_knots <- knots
      modeldata$.metainfo$B <- B
      modeldata$bs_n_basis <- ncol(B)
      modeldata$bs_dists <- c(
        rep( # left boundary basis functions
          knots$interior[1]-knots$boundary[1], spline_degree - 1
        ),
        diff( # interior basis functions
          knots$interior
        ),
        rep( # right boundary basis function
          knots$boundary[2]-knots$interior[length(knots$interior)], 1
        )
      )

      B_sparse <- suppressMessages(rstan::extract_sparse_parts(B))
      modeldata$bs_n_w <- length(B_sparse$w)
      modeldata$bs_w <- B_sparse$w
      modeldata$bs_v <- B_sparse$v
      modeldata$bs_u <- B_sparse$u

      modeldata$.init$bs_coeff_noise_raw <- rep(0, sum(modeldata$bs_dists))
    },
    required = c(
      ".metainfo$length_I",
      ".metainfo$length_seeding",
      ".metainfo$partial_window",
      ".metainfo$forecast_horizon"
      ),
    modeldata = modeldata
  )

  modeldata$inf_ar_sd_prior <- set_prior(
    "inf_ar_sd",
    "truncated normal",
    mu = inf_sd_prior_mu,
    sigma = inf_sd_prior_sigma
  )
  modeldata$.init$inf_ar_sd <- init_from_location_scale_prior(
    modeldata$inf_ar_sd_prior
  )

  modeldata$inf_smooth <- inf_smooth
  modeldata$inf_trend_smooth <- inf_trend_smooth
  modeldata$inf_trend_dampen <- inf_trend_dampen

  modeldata$R_w <- R_window

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_approx = c()
  )

  return(modeldata)
}


#' Add change point model for the variability of Rt
#'
#' @param length_R Length of the modeled Rt time series
#' @param changepoint_dist How many days should the change points be apart? If
#'   zero, no change points are modeled (fixed R variability).
#' @inheritParams template_model_helpers
#'
#' @details The change points are placed going backwards from the end of the
#'   time series, i.e. on the last day, then `changepoint_dist` days before
#'   that, and so on. The distance between the first and second change point at
#'   the start of the time series (partly falls into the seeding phase) varies
#'   between half and double the `changepoint_dist`.
#'
#' @inherit modeldata_init return
#' @keywords internal
add_R_variability <- function(length_R, h, length_seeding, length_partial,
                              length_partial_gen,
                              changepoint_dist, modeldata) {
  if (changepoint_dist == 0) {
    # no changepoints
    B <- matrix(rep(1,length_R+h), ncol = 1)
  } else if (length_R - length_partial - length_partial_gen - changepoint_dist <= length_seeding) {
    B <- matrix(rep(1,length_R+h), ncol = 1)
  } else {
    # we here suppress extrapolation warnings if h > changepoint_dist, as we fix
    # spline bases in the next step
    B <- suppressWarnings(splines::bs(
      1:(length_R+h),
      knots = rev(c(seq(
        length_R - length_partial - length_partial_gen,
        length_seeding + changepoint_dist/2,
        by = -changepoint_dist),
        length_seeding
      )),
      degree = 1,
      intercept = TRUE,
      Boundary.knots = c(1, length_R),
    ))
    B <- matrix(pmax(0, pmin(B, 1)), ncol = ncol(B)) # fix extrapolated bases
    if (all(B[,ncol(B)]==0)) {
      B <- B[,-ncol(B)]
    }
  }
  modeldata$R_vari_n_basis <- ncol(B)
  B_sparse <- suppressMessages(rstan::extract_sparse_parts(B))
  modeldata$R_vari_n_w <- length(B_sparse$w)
  modeldata$R_vari_w <- B_sparse$w
  modeldata$R_vari_v <- B_sparse$v
  modeldata$R_vari_u <- B_sparse$u
  return(modeldata)
}

#' Estimate constant seeding infections
#'
#' @description This option estimates a constant number of infections at the
#'   start of the modeled time period.
#'
#' @inheritParams seeding_estimate_rw
#'
#' @details The seeding phase has the length of the maximum generation time
#'   (during this time, the renewal model cannot be applied). It is here assumed
#'   that the expected number of new infections stays constant over this time
#'   period. This assumption can however be violated: While traditionally,
#'   seeding refers to the first few (potentially imported) infections of an
#'   epidemic, depending on what time period the model is fitted to, it may also
#'   cover a different phase with stronger growth dynamics. Thus, if your data
#'   starts in the middle of an epidemic wave, it is recommended to use
#'   [seeding_estimate_rw()] instead of [seeding_estimate_constant()].
#'
#' @details If `intercept_prior_q5` or `intercept_prior_q95` are not specified
#'   by the user, `EpiSewer` will compute a rough median empirical estimate of
#'   the number of cases using the supplied wastewater measurements and shedding
#'   assumptions, and then infer the missing quantiles based on this. If none of
#'   the quantiles are provided, they are set to be roughly 1/10 and 10 times
#'   the empirical median estimate. We note that this is a violation of Bayesian
#'   principles (data must not be used to inform priors) - but a neglectable
#'   one, since it only ensures that the seeding is modeled on the right order
#'   of magnitude and does not have relevant impacts on later Rt estimates.
#'
#' @details The priors of this component have the following functional form:
#' - initial number of infections (log scale): `Normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
seeding_estimate_constant <- function(
    intercept_prior_q5 = NULL,
    intercept_prior_q95 = NULL,
    extend = TRUE,
    modeldata = modeldata_init()) {

  modeldata$seeding_model <- 0

  modeldata <- add_seeding_intercept_prior(
    intercept_prior_q5 = intercept_prior_q5,
    intercept_prior_q95 = intercept_prior_q95,
    calling_f_name = "seeding_estimate_constant",
    modeldata = modeldata
  )

  modeldata$iota_log_seed_sd_prior <- numeric(0)
  modeldata$.init$iota_log_seed_sd <- numeric(0)
  modeldata$.init$iota_log_ar_noise <- numeric(0)

  modeldata$iota_log_seed_trend_reg <- tbe(
    rep(0, modeldata$.metainfo$length_seeding),
    ".metainfo$length_seeding"
  )

  modeldata$.metainfo$extend_seeding <- extend

  modeldata$.str$infections[["seeding"]] <- list(
    seeding_estimate_fixed = c()
  )

  return(modeldata)
}

#' Estimate seeding infections using a random walk model
#'
#' @description This option estimates initial infections at the start of the
#'   modeled time period, when the renewal model cannot be applied yet. It uses
#'   a geometric random walk to model these seeding infections.
#'
#' @param intercept_prior_q5 Prior (5% quantile) on the initial number of
#'   infections. Can be interpreted as an approximate lower bound. If NULL
#'   (default), this is computed from a crude empirical estimate of the number
#'   of cases (see details).
#' @param intercept_prior_q95 Prior (95% quantile) on the initial number of
#'   infections. Can be interpreted as an approximate upper bound. If NULL
#'   (default), this is computed from a crude empirical estimate of the number
#'   of cases (see details).
#' @param rel_change_prior_mu Prior (mean) on the relative change rate of the
#'   geometric random walk during the seeding phase. The default value (0.05)
#'   assumes that daily changes are +-5% on expectation and likely less than
#'   +-10% per day.
#' @param rel_change_prior_sigma Prior (standard deviation) on the relative
#'   change rate of the geometric random walk during the seeding phase. This
#'   expresses your uncertainty about the change rate. The default value (0.025)
#'   assumes that the daily change rate could be 5% points higher or lower than
#'   your prior mean. For example, if `rel_change_prior_mu = 0.05` and
#'   `rel_change_prior_sigma = 0.025`, this means you expect the daily change
#'   rate to be between 0 (0%) and 0.1 (10%).
#' @param extend Should the seeding phase be extended when concentrations are
#'   very low at the start of the measurement time series? If `TRUE`, then the
#'   seeding phase will be extended to the first date with three consecutive
#'   detects (i.e. non-zero measurements). The reproduction number will only be
#'   modeled from that date onward. This option often makes sense, as infection
#'   numbers are typically very low during a period with many non-detects, which
#'   can lead to sampling problems when estimating Rt. If you nevertheless want
#'   Rt estimates also for this period, you can use `extend = FALSE`. Note
#'   though that estimated reproduction numbers are not necessarily meaningful
#'   during periods with very low infection numbers, as transmission dynamics
#'   may be dominated by chance events and importations.
#'
#' @details The seeding phase has the length of the maximum generation time
#'   (during this time, the renewal model cannot be applied). Traditionally,
#'   seeding refers to the first few (potentially imported) infections of an
#'   epidemic, but depending on what time period the model is fitted to, this
#'   may also cover a different phase with stronger growth dynamics.
#'
#' @details If `intercept_prior_q5` or `intercept_prior_q95` are not specified
#'   by the user, `EpiSewer` will compute a rough median empirical estimate of
#'   the number of cases using the supplied wastewater measurements and shedding
#'   assumptions, and then infer the missing quantiles based on this. If none of
#'   the quantiles are provided, they are set to be roughly 1/10 and 10 times
#'   the empirical median estimate. We note that this is a violation of Bayesian
#'   principles (data must not be used to inform priors) - but a neglectable
#'   one, since it only ensures that the seeding is modeled on the right order
#'   of magnitude and does not have relevant impacts on later Rt estimates.
#'
#' @details The priors of this component have the following functional form:
#' - intercept of the random walk (log scale): `Normal`
#' - standard deviation of the random walk (log scale): `Truncated normal`
#'   The priors for these parameters are determined based on the user-supplied
#'   arguments, using appropriate transformations and the two-sigma-rule of
#'   thumb.
#'
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
seeding_estimate_rw <- function(
    intercept_prior_q5 = NULL,
    intercept_prior_q95 = NULL,
    rel_change_prior_mu = 0.05,
    rel_change_prior_sigma = 0.025,
    extend = TRUE,
    modeldata = modeldata_init()) {

  modeldata$seeding_model <- 1

  modeldata <- add_seeding_intercept_prior(
    intercept_prior_q5 = intercept_prior_q5,
    intercept_prior_q95 = intercept_prior_q95,
    calling_f_name = "seeding_estimate_rw",
    modeldata = modeldata
    )

  modeldata$iota_log_seed_sd_prior <- set_prior(
    "iota_log_seed_sd",
    "truncated normal",
    mu = rel_change_prior_mu,
    sigma = rel_change_prior_sigma
  )

  modeldata$.init$iota_log_seed_sd <- 1
  modeldata$.init$iota_log_ar_noise <- tbe(
    rep(0, modeldata$.metainfo$length_seeding - 1),
    ".metainfo$length_seeding"
  )

  # compute regression vector for estimating log-linear trend of seeding phase
  modeldata$iota_log_seed_trend_reg <- tbe(
    get_regression_linear_trend(
      1:modeldata$.metainfo$length_seeding,
      weights = c(rep(0,modeldata$se), rev(modeldata$generation_dist))
      ),
    c(".metainfo$length_seeding", "generation_dist")
  )

  modeldata$.metainfo$extend_seeding <- extend

  modeldata$.str$infections[["seeding"]] <- list(
    seeding_estimate_rw = c()
  )

  return(modeldata)
}

add_seeding_intercept_prior <- function(
    intercept_prior_q5, intercept_prior_q95,
    calling_f_name,
    modeldata = modeldata_init()) {

  help_seeding_f <- cli_help(calling_f_name)

  if (!is.null(intercept_prior_q5) && intercept_prior_q5 < 1) {
    cli::cli_warn(paste0(
      "Warning from ",
      help_seeding_f, ": ",
      "The 5% quantile (`intercept_prior_q5`) is very low ",
      "(less than 1 infection)."
    ))
  }

  if (!is.null(intercept_prior_q95) && intercept_prior_q95 < 2) {
    cli::cli_warn(paste0(
      "Warning from ",
      help_seeding_f, ": ",
      "The 95% quantile (`intercept_prior_q95`) is very low ",
      "(less than 2 infections)."
    ))
  }

  if (!is.null(intercept_prior_q5) && !is.null(intercept_prior_q95)) {
    modeldata$iota_log_seed_intercept_prior <- set_prior_normal_log(
      "iota_log_seed_intercept",
      unit_q5 = intercept_prior_q5,
      unit_q95 = intercept_prior_q95
    )
    modeldata$.init$iota_log_seed_intercept <- init_from_location_scale_prior(
      modeldata$iota_log_seed_intercept_prior
    )
  } else {
    modeldata <- tbc(
      "seeding_prior",
      {
        cases_crude <- modeldata$.metainfo$initial_cases_crude
        if (!is.null(intercept_prior_q5) && intercept_prior_q5 > cases_crude) {
          cli::cli_warn(paste0(
            "Warning from ",
            help_seeding_f, ": ",
            "You specified a lower quantile ",
            paste0("(`intercept_prior_q5=", intercept_prior_q5, "`) "),
            "that seems be quite large based on the data provided ",
            "(greater than the crude median estimate of ",
            round(cases_crude), " cases). ",
            "The median is now assumed to be 2 times your 5% quantile. ",
            "Please consider specifying a smaller 5% quantile or adjust other ",
            "assumptions (e.g. ",
            cli_help("load_per_case_assume", "load_per_case"),
            " if the crude median estimate seems wrong."
          ))
          cases_crude <- intercept_prior_q5 * 2
        }
        if (!is.null(intercept_prior_q95) && intercept_prior_q95 < cases_crude) {
          cli::cli_warn(paste0(
            "Warning from ",
            help_seeding_f, ": ",
            "You specified an upper quantile ",
            paste0("(`intercept_prior_q95=", intercept_prior_q95, "`) "),
            "that seems be quite low based on the data provided ",
            "(lower than the crude median estimate of ",
            round(cases_crude), " cases). ",
            "The median is now assumed to be half of your 95% quantile. ",
            "Please consider specifying a higher 95% quantile or adjust other ",
            "assumptions (e.g. ",
            cli_help("load_per_case_assume", "load_per_case"),
            " if the crude median estimate seems wrong."
          ))
          cases_crude <- intercept_prior_q95 / 2
        }
        modeldata$iota_log_seed_intercept_prior <- set_prior_normal_log(
          "iota_log_seed_intercept",
          unit_median = cases_crude,
          unit_q5 = intercept_prior_q5,
          unit_q95 = intercept_prior_q95,
          unit_factor = 10,
          paste_dist = ", (mu based on crude empirical estimate of cases)"
        )
        modeldata$.init$iota_log_seed_intercept <- init_from_location_scale_prior(
          modeldata$iota_log_seed_intercept_prior
        )
      },
      required = ".metainfo$initial_cases_crude",
      modeldata = modeldata
    )
  }
  return(modeldata)
}

#' Do not model infection noise
#'
#' This option does not model noise in the infection process and instead
#' implements a deterministic renewal model.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {infection noise models}
infection_noise_none <- function(modeldata = modeldata_init()) {
  modeldata$I_sample <- FALSE
  modeldata$I_overdispersion <- FALSE
  modeldata$I_xi_prior <- numeric(0)
  modeldata$.init$I_xi <- numeric(0)
  modeldata$.init$I <- numeric(0)
  modeldata$.init$I_log <- numeric(0)

  modeldata$.str$infections[["infection_noise"]] <- list(
    infection_noise_none = c()
  )

  return(modeldata)
}

#' Estimate infection noise
#'
#' @description This option estimates noise in the infection process, i.e.
#'   implements a stochastic renewal model. This allows for variation in the
#'   number of new infections generated at each time step, which can often speed
#'   up model fitting.
#'
#' @param overdispersion If `TRUE` (default), new infections are modeled as
#'   Negative Binomial distributed. If `FALSE`, new infections are modeled as
#'   Poisson distributed.
#' @param overdispersion_prior_mu Prior (mean) on the overdispersion parameter
#'   of the Negative Binomial. The default of 0.1 corresponds to 10%
#'   overdispersion. It is also the limit of the coefficient of variation (CV)
#'   of infections as the infection incidence becomes large.
#' @param overdispersion_prior_sigma Prior (standard deviation) on the
#'   overdispersion parameter of the Negative Binomial. If this is set to zero
#'   (default), the overdispersion will be fixed to the prior mean and not
#'   estimated.
#'
#' @details The level of overdispersion is often unidentifiable from a single
#'   time series of measurements. This is why the overdispersion is fixed by
#'   default.
#'
#' @details For complicated reasons, MCMC sampling of high infection numbers is
#'   faster with a certain degree of overdispersion. Thus, if modeling large
#'   waves, assuming some overdispersion can also make sense for computational
#'   reasons. The effects on the estimated transmission dynamics are often
#'   minimal.
#'
#' @details The priors of this component have the following functional form:
#' - overdispersion parameter of the Negative Binomial: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {infection noise models}
infection_noise_estimate <-
  function(overdispersion = TRUE,
           overdispersion_prior_mu = 0.1,
           overdispersion_prior_sigma = 0,
           modeldata = modeldata_init()) {

    if (overdispersion_prior_mu == 0 && overdispersion_prior_sigma == 0) {
      modeldata$I_overdispersion <- FALSE
    } else {
      modeldata$I_overdispersion <- overdispersion
    }

    if (modeldata$I_overdispersion) {
      modeldata$I_xi_prior <- set_prior(
        "I_xi", "normal",
        mu = overdispersion_prior_mu, sigma = overdispersion_prior_sigma
      )
      modeldata$.init$I_xi <- init_from_location_scale_prior(modeldata$I_xi_prior)
    } else {
      modeldata$I_xi_prior <- numeric(0)
      modeldata$.init$I_xi <- numeric(0)
    }

    modeldata$I_sample <- TRUE
    modeldata$.init$I <- tbe(
      rep(
        modeldata$.metainfo$initial_cases_crude,
        modeldata$.metainfo$length_I
      ),
      c(".metainfo$initial_cases_crude", ".metainfo$length_I")
    )
    modeldata$.init$I_log <- tbe(
      rep(
        log(modeldata$.metainfo$initial_cases_crude),
        modeldata$.metainfo$length_I
      ),
      c(".metainfo$initial_cases_crude", ".metainfo$length_I")
    )
    modeldata$.init$I_raw <- tbe(
      rep(
        0,
        modeldata$.metainfo$length_I
      ),
      c(".metainfo$length_I")
    )

    modeldata$.str$infections[["infection_noise"]] <- list(
      infection_noise_estimate = c(overdispersion = overdispersion)
    )

    return(modeldata)
  }

add_link_function <- function(link, R_max, modeldata) {
  # link function and hyperparameters
  if (link == "inv_softplus") {
    # first element: 0 = inv_softplus
    # second element: sharpness of softplus function
    # the last two entries are just placeholders for hyperparameters of
    # other link functions
    modeldata$R_link = c(0, 4, 0, 0)
  } else if (link == "scaled_logit") {
    # logistic function for scaled logit link
    if (R_max < 3) {
      cli::cli_abort("Please chose a maximum R value greater or equal to 3.")
    }
    modeldata$.metainfo$R_max <- R_max
    a_k <- logistic_find_a_k(R_max)
    # first element: 1 = scaled_logit
    # hyperparameters are chosen such that R=1 at x=1, and is roughly
    # linear with slope 1 around R=1
    modeldata$R_link <- c(1, R_max, a_k$a, a_k$k)
  } else {
    cli::cli_abort(c(
      paste0("'", link, "' is not a valid link function. ",
             "Please use one of:"), "inv_softplus", "scaled_logit"))
  }
  return(modeldata)
}

get_selection_vector <- function(from, to, n) {
  selection <- rep(0, n)
  selection[max(1,min(n,from)):max(1,min(n,to))] <- 1
  if (from < 1) {
    selection[1] <- selection[1] + (min(1, to)-from)
  }
  if (to > n) {
    selection[n] <- selection[n] + (to-max(n, from))
  }
  return(selection)
}
