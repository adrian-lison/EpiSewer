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
#'   smoothing prior (Gaussian process by default). A variety of smoothing
#'   priors is supported. Modeling options:
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
    R = R_estimate_gp(),
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
#'  damping through an innovations state space model with a level, trend, and
#'  damping component.
#'
#'@param R_start_prior_mu Prior (mean) on the initial value of Rt (level).
#'@param R_start_prior_sigma Prior (standard deviation) on the initial value of
#'  Rt (level).
#'@param trend_prior_mu Prior (mean) on the initial trend of Rt.
#'@param trend_prior_sigma Prior (standard deviation) on the initial trend of
#'  Rt.
#'@param sd_base_prior_mu Prior (mu) on the baseline standard deviation of the
#'  innovations. Please note that for consistency, the overall standard
#'  deviation of innovations will always be the baseline plus an additive
#'  component from `sd_change_prior` - even if no changepoints are modeled (see
#'  below).
#'@param sd_base_prior_sd Prior (standard deviation) on the baseline standard
#'  deviation of the innovations. See `sd_base_prior_mu` for details.
#'@param sd_change_distance Distance between changepoints used to model
#'  additional variation in Rt. The default change point distance is 4 weeks.
#'  Very short changepoint distances must be chosen with care, as they can make
#'  the Rt time series too flexible. If set to zero, no change points are
#'  modeled.
#'@param sd_change_prior_shape Exponential-Gamma prior (shape) on standard
#'  deviation additional to baseline. This prior describes the distribution of
#'  the standard deviation of Rt over time. EpiSewer will estimate a baseline
#'  standard deviation (see `sd_base_prior_mu`), and model additional variation
#'  on top of the baseline using a changepoint model. Please see the details for
#'  more explanation.
#'@param sd_change_prior_scale Exponential-Gamma prior (scale) on standard
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
#'@param dampen_prior_mu Prior (mean) on the damping parameter. Must be
#'  between 0 and 1.
#'@param dampen_prior_sigma Prior (standard deviation) on the damping
#'  parameter. If this is set to zero, the damping parameter will be fixed to
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
#'  level, a trend, and a damping component.
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
#' - The damping determines how long a previous trend continues into the
#'  future before it levels of to a stationary time series. The strength of
#'  damping is controlled by a damping parameter (often called phi). Note
#'  that *smaller* values of `phi` indicate *stronger* damping. In particular,
#'  `phi = 1` means no damping. Values below `phi = 0.8` are seldom in
#'  practice as the damping becomes very strong.
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
#'  shape and scale of this Gamma distribution. The distribution has a strong
#'  peak towards zero and a long tail. This regularizes the estimated deviations
#'  from the baseline standard deviation - most deviations are small, but during
#'  special time periods, the deviation might also be larger.
#'
#'@details The priors of this component have the following functional form:
#' - initial level of Rt: `Normal`
#' - initial trend of Rt: `Normal`
#' - baseline standard deviation of innovations: `Half-normal`
#' - additional standard deviation at changepoints: `Exponential-Gamma`
#' - smoothing parameter: `Beta`
#' - trend smoothing parameter: `Beta`
#' - damping parameter: `Beta`
#'
#'@inheritParams template_model_helpers
#'@inherit modeldata_init return
#'@export
#'@family {Rt models}
R_estimate_ets <- function(
    R_start_prior_mu = 1,
    R_start_prior_sigma = 0.8,
    trend_prior_mu = 0,
    trend_prior_sigma = 0.1,
    sd_base_prior_mu = 0,
    sd_base_prior_sd = 0.025,
    sd_change_prior_shape = 0.5,
    sd_change_prior_scale = 1e-4,
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

  modeldata <- configure_R_model(
    name_approach = "ets",
    model_id = 0,
    use_ets = TRUE,
    modeldata = modeldata
  )

  modeldata$R_intercept_prior <- set_prior(
    "R_intercept", "normal",
    mu = R_start_prior_mu, sigma = R_start_prior_sigma
  )

  modeldata$ets_trend_start_prior <- set_prior(
    "ets_trend_start", "normal",
    mu = trend_prior_mu, sigma = trend_prior_sigma
  )

  modeldata$R_sd_baseline_prior <- set_prior("R_sd_baseline",
     "normal",
     mu = sd_base_prior_mu,
     sigma = sd_base_prior_sd
  )
  modeldata$R_sd_change_prior <- set_prior("R_sd_change",
    "lomax",
    shape = sd_change_prior_shape,
    scale = sd_change_prior_scale
  )

  modeldata$.init$R_intercept <-
    modeldata$R_intercept_prior$R_intercept_prior[1]
  modeldata$.init$ets_trend_start <- 1e-4

  modeldata <- tbc(
    "R_ets_noise",
    {
      modeldata$ets_length <- modeldata$.metainfo$length_R_modeled
      modeldata$.init$ets_noise <- rep(0, modeldata$.metainfo$length_R_modeled - 1)
      modeldata <- add_R_variability(
        length_R = modeldata$.metainfo$length_R_modeled,
        length_seeding = modeldata$.metainfo$length_seeding,
        partial_window = modeldata$.metainfo$partial_window,
        partial_generation = modeldata$.metainfo$partial_generation,
        h = modeldata$.metainfo$forecast_horizon,
        changepoint_dist = sd_change_distance,
        modeldata = modeldata
      )
      modeldata$.init$R_sd_baseline <- 1e-2
      if (modeldata$R_vari_ncol > 2) {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_ncol - 2)
      } else {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_ncol)
      }
    },
    required = c(
      ".metainfo$length_R_modeled",
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

  # dummies
  modeldata <- add_dummies_basis_splines(modeldata)
  modeldata <- add_dummies_basis_splines2(modeldata)
  modeldata <- add_dummies_R_vari_selection(modeldata)
  modeldata <- add_dummies_soft_changepoints(modeldata)
  modeldata <- add_dummies_smooth_derivative(modeldata)
  modeldata <- add_dummies_gp(modeldata)

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
#'@param R_start_prior_mu Prior (mean) on the initial value of Rt.
#'@param R_start_prior_sigma Prior (standard deviation) on the initial value of
#'  Rt.
#'@param sd_base_prior_mu Prior (mean) on the baseline standard deviation of the
#'  innovations. Please note that for consistency, the overall standard
#'  deviation of innovations will always be the baseline plus an additive
#'  component from `sd_change_prior` even if no changepoints are modeled (see
#'  below).
#'@param sd_base_prior_sd Prior (standard deviation) on the baseline standard
#'  deviation of the innovations. See `sd_base_prior_mu` for details.
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
#'@param sd_change_prior_scale Exponential-Gamma prior (scale) on standard
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
#'  shape and scale of this Gamma distribution. The distribution has a strong
#'  peak towards zero and a long tail. This regularizes the estimated deviations
#'  from the baseline standard deviation - most deviations are small, but during
#'  special time periods, the deviation might also be larger.
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
    R_start_prior_mu = 1,
    R_start_prior_sigma = 0.8,
    sd_base_prior_mu = 0,
    sd_base_prior_sd = 0.025,
    sd_change_prior_shape = 0.5,
    sd_change_prior_scale = 1e-4,
    sd_change_distance = 7*26,
    link = "inv_softplus",
    R_max = 6,
    differenced = FALSE,
    noncentered = TRUE,
    modeldata = modeldata_init()) {
  modeldata <- R_estimate_ets(
    R_start_prior_mu = R_start_prior_mu,
    R_start_prior_sigma = R_start_prior_sigma,
    sd_base_prior_mu = sd_base_prior_mu,
    sd_base_prior_sd = sd_base_prior_sd,
    sd_change_prior_shape = sd_change_prior_shape,
    sd_change_prior_scale = sd_change_prior_scale,
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
#'@description This option estimates the effective reproduction number Rt using
#'  penalized cubic basis splines.
#'
#'@param knot_distance_global Distance between spline breakpoints (knots) for
#'  the *global* spline in days. `EpiSewer` uses an ensemble of two splines to
#'  model Rt. The global spline models larger, long-term changes. Default is
#'  4*7, i.e. one knot every four week.
#'@param knot_distance_local Distance between spline breakpoints (knots) for the
#'  *local* spline in days. `EpiSewer` uses an ensemble of two splines to model
#'  Rt. The local spline models smaller, short-term changes. Default is 7, i.e.
#'  one knot each week.
#'@param R_start_prior_mu Prior (mean) on the initial reproduction number
#'  (intercept).
#'@param R_start_prior_sigma Prior (standard deviation) on the initial
#'  reproduction number (intercept).
#'@param R_sd_local_prior_mu Prior (mean) for the variation of the local (i.e.
#'  short-term) spline. This controls the standard deviation of a random walk
#'  over the coefficients of the local spline. The prior refers to the *daily*
#'  standard deviation and is thus independent of the knot distance.
#'@param R_sd_local_prior_sd Prior (standard deviation) for the variation of the
#'  local (i.e. short-term) spline. See `R_sd_local_prior_mu` for details.
#'@param R_sd_global_prior_shape Exponential-Gamma prior (shape) for the
#'  variation of the global (i.e. long-term) spline. This controls the standard
#'  deviation of a random walk over the coefficients of the global spline. The
#'  prior refers to the *daily* standard deviation and is thus independent of
#'  the knot distance. The exponential-Gamma prior is sparse, i.e. it has a
#'  strong peak towards zero and a long tail. Smaller shape parameters will lead
#'  to more sparseness, i.e. a longer tail. Note that when adjusting the shape,
#'  you will likely also have to adjust the scale. The variation of the global
#'  splines follows a change point model to allow for adaptive changes, see
#'  details for more explanation.
#'@param R_sd_global_prior_scale Exponential-Gamma prior (scale) for the
#'  variation of the global (i.e. long-term) spline. Larger scales will lead to
#'  more variability. See `R_sd_global_prior_shape` and the details for more
#'  explanation about this prior.
#'@param R_sd_global_change_distance Distance between changepoints used to model
#'  global variation in Rt. `EpiSewer` uses an adaptive model for the variation
#'  of the global spline, to model both time periods with stable and with
#'  volatile transmission dynamics. The default change point distance is equal
#'  to `knot_distance_global`. Making this smaller than the global knot distance
#'  will not have a large effect, however making this longer reduces the
#'  adaptability of the Rt variation. If set to zero, no change points are
#'  modeled, meaning zero adaptability.
#'@param link Link function. Currently supported are `inv_softplus` (default)
#'  and `scaled_logit`. Both of these links are configured to behave
#'  approximately like the identity function around R=1, but become increasingly
#'  non-linear below (and in the case of `scaled_logit` also above) R=1.
#'@param R_max If `link=scaled_logit` is used, a maximum reproduction number
#'  must be assumed. This should be higher than any realistic R value for the
#'  modeled pathogen. Default is 6.
#'
#'@details `EpiSewer` uses a combination of two splines (global, for
#'  larger/long-term changes and local, for smaller/short-term changes), with a
#'  random walk on the coefficients of each spline. This allows a highly
#'  adaptive yet regularized Rt model.
#'   - The prior `R_start_prior` should reflect
#'  your expectation of Rt at the beginning of the time series. If your time
#'  series starts in the midst of an outbreak, you might want to use a prior
#'  with mean larger than 1 for the intercept.
#'   - The prior `R_sd_local_prior` on the standard deviation of the local
#'  spline coefficients should be interpreted in terms of daily additive changes
#'  (this is accurate around Rt=1, and becomes less accurate as Rt approaches 0
#'  or its upper bound as defined by the `link` function). For example, a
#'  baseline half-normal prior with sd=0.05 allows a daily standard deviation
#'  between 0 and 0.1. A daily standard deviation of 0.1 in turn roughly allows
#'  the spline coefficients to change by ±0.2 (using the 2 sigma rule) each day.
#'  The daily standard deviation is summed up over the days between two knots to
#'  get the actual standard deviation of the coefficients. This way, the prior
#'  is independent of the chosen `knot_distance`. For example, if
#'  `knot_distance` is 7 days, and a constant daily standard deviation of 0.1 is
#'  estimated, the coefficients of two adjacent splines can differ by up to
#'  `0.2*sqrt(knot_distance)`, i.e. ±0.5. Note however that this does not
#'  directly translate into a change of Rt by ±0.5, as Rt is always the weighted
#'  sum of several basis functions at any given point. It will therefore change
#'  more gradually, depending on the distances between knots.
#'  - The prior `R_sd_global_prior` on the standard deviation of the global
#'  spline coefficients is interpreted in the same way as the local spline, i.e.
#'  it describes daily additive changes. However, we us a *sparse* prior for the
#'  global variation, combined with a changepoint model. This yields a global
#'  spline that expects small changes most of the time but allows for occasional
#'  large changes. For example, during the height of an epidemic wave,
#'  countermeasures may lead to much faster changes in Rt than observable at
#'  other times. These differences in variability are accounted for using change
#'  points placed at regular intervals, with the daily global standard deviation
#'  evolving linearly between the change points. The values at the changepoints
#'  are modeled as independently distributed and follow the Exponential-Gamma
#'  (EG) prior distribution defined by `R_sd_global_prior`. The EG distribution
#'  is also known as Lomax distribution and corresponds to an exponential
#'  distribution with a Gamma distributed rate parameter.
#'
#'@details The smoothness of the Rt estimates is influenced by a combination of
#'  the global and local knot distances and by the priors on the global and
#'  local standard deviation of the random walk on spline coefficients. Placing
#'  knots further apart increases the smoothness of Rt estimates and can speed
#'  up model fitting. The Rt time series remains surprisingly flexible even at
#'  larger knot distances, but placing knots too far apart can lead to
#'  inaccurate estimates. Note that the priors for the variation of the global
#'  and local random walk also influence the uncertainty of Rt estimates towards
#'  the present / date of estimation, when limited data signal is available.
#'  Absent sufficient data signal, Rt estimates will tend to stay at the current
#'  level (which corresponds to assuming unchanged transmission dynamics).
#'
#'@details The priors of this component have the following functional form:
#' - initial Rt (intercept):
#'  `Normal`
#' - daily standard deviation of the random walk over local spline coefficients:
#'  `Half-normal`
#'- daily standard deviation of the random walk over global spline coefficients:
#'  `Exponential-Gamma`
#'
#'@inheritParams template_model_helpers
#'@inherit modeldata_init return
#'@export
#'@family {Rt models}
R_estimate_splines <- function(
    knot_distance_global = 4*7,
    knot_distance_local = 7,
    R_start_prior_mu = 1,
    R_start_prior_sigma = 0.8,
    R_sd_local_prior_mu = 0,
    R_sd_local_prior_sd = 0.05,
    R_sd_global_prior_shape = 1,
    R_sd_global_prior_scale = 1e-2,
    R_sd_global_change_distance = knot_distance_global,
    link = "inv_softplus",
    R_max = 6,
    modeldata = modeldata_init()) {

  modeldata <- configure_R_model(
    name_approach = "splines",
    model_id = 1,
    use_bs = TRUE,
    use_bs2 = TRUE,
    modeldata = modeldata
  )

  spline_degree <- 3 # fixed to cubic splines

  modeldata <- tbc(
    "spline_definition",
    {
      modeldata$.init$R_intercept <- 1

      # Global spline model for Rt
      knots_global <- place_knots(
        length_R = modeldata$.metainfo$length_R_modeled,
        forecast_horizon = modeldata$.metainfo$forecast_horizon,
        knot_distance = knot_distance_global,
        partial_window = modeldata$.metainfo$partial_window,
        partial_generation = modeldata$.metainfo$partial_generation,
        fix_forecast = TRUE
        )
      modeldata$.metainfo$R_knots_global <- knots_global
      modeldata <- use_basis_splines(
        spline_length = with(modeldata$.metainfo, length_R_modeled + forecast_horizon),
        knots = knots_global,
        degree = spline_degree,
        modeldata = modeldata
        )

      # Local spline model for Rt
      knots_local <- place_knots(
        length_R = modeldata$.metainfo$length_R_modeled,
        forecast_horizon = modeldata$.metainfo$forecast_horizon,
        knot_distance = knot_distance_local,
        partial_window = modeldata$.metainfo$partial_window,
        partial_generation = modeldata$.metainfo$partial_generation,
        fix_forecast = FALSE
      )
      modeldata$.metainfo$R_knots_local <- knots_local
      modeldata <- use_basis_splines2(
        spline_length = with(modeldata$.metainfo, length_R_modeled + forecast_horizon),
        knots = knots_local,
        degree = spline_degree,
        modeldata = modeldata
      )

      # Changepoint model for variability of Rt
      modeldata <- add_R_variability(
        length_R = modeldata$.metainfo$length_R_modeled,
        length_seeding = modeldata$.metainfo$length_seeding,
        partial_window = modeldata$.metainfo$partial_window,
        partial_generation = modeldata$.metainfo$partial_generation,
        h = modeldata$.metainfo$forecast_horizon,
        changepoint_dist = R_sd_global_change_distance,
        modeldata = modeldata
        )
      modeldata$.init$R_sd_baseline <- 1e-2
      if (modeldata$R_vari_ncol > 2) {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_ncol - 2)
      } else {
        modeldata$.init$R_sd_changepoints <- rep(1e-2, modeldata$R_vari_ncol)
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
        n = modeldata$bs_length
      ))
      modeldata <- add_sparse_matrix(
        R_vari_selection, "R_vari_sel", modeldata
        )

      all_positions_local <- c(
        rev(knots_local$interior[1] - seq(
          0,
          by = (knots_local$interior[1] - knots_local$boundary[1]),
          length.out = spline_degree
        )), knots_local$interior[-1], knots_local$boundary[2])
      # this is currently not used because we have
      # no changepoint model for the local spline
      R_vari_selection_local <- t(mapply(
        get_selection_vector,
        from = all_positions_local[-length(all_positions_local)] + 1,
        to = all_positions_local[-1],
        n = modeldata$bs2_length
      ))
      modeldata <- add_sparse_matrix(
        R_vari_selection_local, "R_vari_sel_local", modeldata
        )

    },
    required = c(
      ".metainfo$length_R_modeled",
      ".metainfo$length_seeding",
      ".metainfo$partial_window",
      ".metainfo$partial_generation",
      ".metainfo$forecast_horizon"
      ),
    modeldata = modeldata
  )

  modeldata$R_intercept_prior <- set_prior("R_intercept",
    "normal",
    mu = R_start_prior_mu,
    sigma = R_start_prior_sigma
  )
  modeldata$R_sd_baseline_prior <- set_prior("R_sd_baseline",
    "normal",
    mu = R_sd_local_prior_mu,
    sigma = R_sd_local_prior_sd
    )
  modeldata$R_sd_change_prior <- set_prior("R_sd_change",
    "lomax",
    shape = R_sd_global_prior_shape,
    scale = R_sd_global_prior_scale
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  # dummies
  modeldata <- add_dummies_exponential_smoothing(modeldata)
  modeldata <- add_dummies_soft_changepoints(modeldata)
  modeldata <- add_dummies_smooth_derivative(modeldata)
  modeldata <- add_dummies_gp(modeldata)

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
#' @param inf_trend_dampen Trend damping parameter (phi) for infections. The
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
        length_R = with(modeldata$.metainfo, length_I - length_seeding + forecast_horizon),
        forecast_horizon = 0,
        knot_distance = knot_distance,
        partial_generation = 0,
        partial_window = with(modeldata$.metainfo, partial_window + forecast_horizon)
      )
      modeldata$.metainfo$R_knots <- knots
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

      modeldata <- use_basis_splines(
        spline_length = with(modeldata$.metainfo, length_I - length_seeding + forecast_horizon),
        knots = knots,
        degree = spline_degree,
        modeldata = modeldata
      )
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

#' Estimate Rt via piecewise constant changepoint model
#'
#' @description This option models the effective reproduction number Rt over
#'   time as a piecewise constant function. The changepoints of the constant
#'   pieces are estimated from the data using an approximate changepoint model.
#'   This option is suitable when Rt follows a step-like trajectory with
#'   occasional large jumps, e.g. due to strong interventions. It is less
#'   suitable for modeling gradual changes in Rt over time.
#' @param changepoint_max_distance Maximum distance between changepoints in
#'   days. This setting guarantees that two consecutive changes in Rt can be
#'   captured if they are `changepoint_max_distance` days or more apart. Faster
#'   changes are not guaranteed to be captured.
#' @param changepoint_min_distance Minimum distance between changepoints in
#'   days. This setting guarantees that two consecutive changes in Rt are at
#'   least `changepoint_min_distance` days apart. This avoids Rt trajectories
#'   that are unrealistically volatile. For example, a minimum distance of 7
#'   implies that you don't expect Rt to make more than one large jump within a
#'   week.
#' @param change_prior_shape Exponential-Gamma (EG) prior (shape) for the
#'   strength of changes between the pieces. At each estimated changepoint, the
#'   change in Rt is sampled from a normal distribution with standard deviation
#'   given by the EG prior. This prior has a strong peak towards zero and a long
#'   tail. In other words, while most changes are expected to be small, the
#'   prior allows for occasional large jumps in Rt. Smaller shape parameters
#'   will lead to more sparseness, i.e. a longer tail. Note that when adjusting
#'   the shape, you will likely also have to adjust the scale. See details for
#'   more advice on choosing a suitable prior.
#' @param change_prior_scale Exponential-Gamma (EG) prior (scale) for the
#'   strength of changes between the pieces. See `change_prior_shape` above for
#'   an explanation. Larger scales will lead to more variability: a doubling of
#'   the scale roughly corresponds to a doubling of all quantiles of the prior.
#' @param change_tolerance Tolerance for "negligible" changes in Rt. Changes
#'   smaller than `change_tolerance` are ignored by `changepoint_min_distance`,
#'   i.e. they can also occur closer to each other. This tolerance gives the
#'   model more flexibility in placing changepoints with large jumps.
#' @param strictness_alpha The concentration parameter of the Dirichlet prior
#'   for the changepoint positions. Choosing smaller values of
#'   `strictness_alpha` will lead to more discrete changepoints. Note that
#'   choosing small values of `strictness_alpha` can impede MCMC sampling.
#'
#' @details The Exponential-Gamma (EG) prior on the strength of changes is
#'   parameterized via the arguments `change_prior_shape` and
#'   `change_prior_scale`. It has a long tail to support large changes while
#'   keeping the variation low most of the time. The default configuration
#'   should work well in most contexts except for really extreme changes in Rt
#'   over a short time window. To check the quantiles of your prior, you can use
#'   the function [qexpgamma()] with corresponding shape and scale parameters.
#'
#' @details If you need to adjust the overall variation, you can adjust the
#'   `change_prior_scale` parameter. A doubling of the scale roughly corresponds
#'   to a doubling of the quantiles. For example, when the 95% quantile is 0.2
#'   for a given scale and you double that scale, the 95% quantile will be at
#'   0.4.
#'
#' @details If you need to support more extreme changes, you can decrease the
#'   `change_prior_shape` parameter, which will emphasize the long-tail behavior
#'   of the prior. Note however that this will substantially increase all
#'   quantiles of the prior, so you will also have to decrease the scale
#'   parameter to achieve a similar level of day-to-day variation.
#'
#' @inheritParams R_estimate_splines
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
R_estimate_piecewise <- function(
    R_start_prior_mu = 1,
    R_start_prior_sigma = 0.8,
    changepoint_max_distance = 14,
    changepoint_min_distance = 7,
    change_prior_shape = 0.5,
    change_prior_scale = 1e-4,
    change_tolerance = 0.05,
    link = "inv_softplus",
    R_max = 6,
    strictness_alpha = 1,
    modeldata = modeldata_init()
    ) {

  modeldata <-  configure_R_model(
    name_approach = "piecewise",
    model_id = 2,
    use_scp = TRUE,
    modeldata = modeldata
  )

  # Settings currently hidden from the user:
  #
  # @parameter sharpness_boltzmann The temperature parameter of the Boltzmann
  # operator, which is used for a smooth approximation of the maximum function
  # that is applied as a fuzzy OR on the changepoint indicators.
  #
  # @parameter strictness_tol_k The k parameter of the logistic function used to
  #   soft-constrain the change tolerance bound. A minimum value of 4 is
  #   recommended, higher values are likely not necessary and can impede
  #   sampling.
  sharpness_boltzmann <- 25
  strictness_tol_k <- 4

  modeldata$R_intercept_prior <- set_prior(
    "R_intercept", "normal",
    mu = R_start_prior_mu, sigma = R_start_prior_sigma
  )
  modeldata$.init$R_intercept <-
    modeldata$R_intercept_prior$R_intercept_prior[1]

  modeldata$R_sd_change_prior <- set_prior("R_sd_change",
    "lomax",
    shape = change_prior_shape,
    scale = change_prior_scale
  )

  modeldata <- tbc(
    "R_piecewise",
    {
      modeldata$.metainfo$R_knots <- knots
      modeldata <- use_soft_changepoints(
        scp_length = modeldata$.metainfo$length_R_modeled,
        last_knot = modeldata$.metainfo$length_R_modeled - changepoint_min_distance,
        distance = changepoint_max_distance,
        min_distance = changepoint_min_distance,
        min_distance_tolerance = change_tolerance,
        strictness_tol_k = strictness_tol_k,
        sharpness_boltzmann = sharpness_boltzmann,
        strictness_alpha = strictness_alpha,
        modeldata
        )
    },
    required = c(
      ".metainfo$length_R_modeled",
      ".metainfo$partial_window",
      ".metainfo$partial_generation"
    ),
    modeldata = modeldata
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  # dummy data
  modeldata <- add_dummies_exponential_smoothing(modeldata)
  modeldata <- add_dummies_basis_splines(modeldata)
  modeldata <- add_dummies_basis_splines2(modeldata)
  modeldata <- add_dummies_R_vari_selection(modeldata)
  modeldata <- add_dummies_R_vari(modeldata)
  modeldata <- add_dummies_smooth_derivative(modeldata)
  modeldata <- add_dummies_gp(modeldata)

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_piecewise = c()
  )

  return(modeldata)
}

#' Estimate Rt via a changepoint spline model
#'
#' @description This option models the effective reproduction number Rt over
#'   time using cubic splines that are regularized via a linear changepoint
#'   model. The Rt trajectory will thus follow a smoothed linear trend, with
#'   time points of trend changes estimated from the data using an approximate
#'   changepoint model. This approach offers high flexibility of the Rt
#'   trajectory while avoiding overfitting on noise.
#' @param changepoint_max_distance Maximum distance (in days) between changes of
#'   the Rt trend. This setting guarantees that two consecutive changes in the
#'   Rt trend can be captured if they are `changepoint_max_distance` days or
#'   more apart. Faster changes are not guaranteed to be captured.
#' @param changepoint_min_distance Minimum distance (in days) between changes of
#'   the Rt trend. This setting guarantees that two consecutive changes in the
#'   Rt trend are at least `changepoint_min_distance` days apart. This avoids Rt
#'   trajectories that are unrealistically volatile. For example, a minimum
#'   distance of 7 implies that you don't expect Rt to significantly change its
#'   trend more than once within a week.
#' @param trend_prior_shape Exponential-Gamma (EG) prior (shape) for the trend
#'   in Rt. At each estimated changepoint, the Rt trend is sampled from a normal
#'   distribution with standard deviation given by the EG prior. This prior has
#'   a strong peak towards zero and a long tail. In other words, while we expect
#'   Rt to remain stable most of the time, this prior also allows for occasional
#'   strong trends. Smaller shape parameters will lead to a longer tail, hence
#'   more extreme trends are supported. Note that when adjusting the shape, you
#'   will likely also have to adjust the scale. See details for more advice on
#'   choosing a suitable prior.
#' @param trend_prior_scale Exponential-Gamma (EG) prior (scale) for the trend
#'   in Rt. See `change_prior_shape` above for an explanation. Larger scales
#'   will lead to more variability: a doubling of the scale roughly corresponds
#'   to a doubling of all quantiles of the prior.
#' @param trend_change_tolerance Tolerance for "negligible" trend changes.
#'   Differences in the trend that are smaller than `change_tolerance` are
#'   ignored by `changepoint_min_distance`, i.e. they can also occur closer to
#'   each other. This tolerance gives the model more flexibility in placing
#'   changepoints with large trend changes.
#' @param strictness_alpha The concentration parameter of the Dirichlet prior
#'   for the changepoint positions. Choosing smaller values of
#'   `strictness_alpha` will lead to more strict changepoints. Note that
#'   choosing small values of `strictness_alpha` can impede MCMC sampling.
#'
#' @details The Exponential-Gamma (EG) prior on the Rt trend is parameterized
#'   via the arguments `change_prior_shape` and `change_prior_scale`. It has a
#'   long tail to support large trends while keeping the Rt variation low most
#'   of the time. The default configuration should work well in most contexts
#'   except for really extreme changes in Rt over a short time window. To check
#'   the quantiles of your prior, you can use the function [qexpgamma()] with
#'   corresponding shape and scale parameters.
#'
#' @details If you need to adjust the overall variation, you can adjust the
#'   `change_prior_scale` parameter. A doubling of the scales roughly
#'   corresponds to a doubling of the quantiles. For example, when the 95%
#'   quantile is 0.2 for a given scale and you double that scale, the 95%
#'   quantile will be at 0.4.
#'
#' @details If you need to support more extreme changes, you can decrease the
#'   `change_prior_shape` parameter, which will emphasize the long-tail behavior
#'   of the prior. Note however that this will substantially increase all
#'   quantiles of the prior, so you will also have to decrease the scale
#'   parameter to achieve a similar level of day-to-day variation.
#'
#' @inheritParams R_estimate_splines
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
R_estimate_changepoint_splines <- function(
    R_start_prior_mu = 1,
    R_start_prior_sigma = 0.8,
    changepoint_max_distance = 3*5,
    changepoint_min_distance = 3*2,
    trend_prior_shape = 5e1,
    trend_prior_scale = 10e-1,
    trend_change_tolerance = 0.01,
    spline_knot_distance = 3,
    link = "inv_softplus",
    R_max = 6,
    strictness_alpha = 0.5,
    modeldata = modeldata_init()
) {

  modeldata <-  configure_R_model(
    name_approach = "changepoint_splines",
    model_id = 3,
    use_bs = TRUE,
    use_scp = TRUE,
    modeldata = modeldata
  )

  # Settings currently hidden from the user:
  #
  # @parameter sharpness_boltzmann The temperature parameter of the Boltzmann
  # operator, which is used for a smooth approximation of the maximum function
  # that is applied as a fuzzy OR on the changepoint indicators.
  #
  # @parameter strictness_tol_k The k parameter of the logistic function used to
  #   soft-constrain the change tolerance bound. A minimum value of 4 is
  #   recommended, higher values are likely not necessary and can impede
  #   sampling.
  sharpness_boltzmann <- 25
  strictness_tol_k <- 4

  modeldata$R_intercept_prior <- set_prior(
    "R_intercept", "normal",
    mu = R_start_prior_mu, sigma = R_start_prior_sigma
  )
  modeldata$.init$R_intercept <-
    modeldata$R_intercept_prior$R_intercept_prior[1]

  modeldata$R_sd_change_prior <- set_prior("R_sd_change",
                                           "lomax",
                                           shape = trend_prior_shape,
                                           scale = trend_prior_scale
  )

  modeldata <- tbc(
    "R_changepoint_splines",
    {
      # Spline model
      spline_knots <- list(
        interior = rev(seq(
          modeldata$.metainfo$length_R_modeled - spline_knot_distance,
          1, by = -spline_knot_distance
        )),
        boundary = c(-3, modeldata$.metainfo$length_R_modeled)
      )
      modeldata$.metainfo$R_knots <- spline_knots
      modeldata <- use_basis_splines(
        spline_length = modeldata$.metainfo$length_R_modeled,
        knots = spline_knots,
        degree = 3,
        modeldata = modeldata
      )

      # Changepoint model
      last_change <- with(modeldata$.metainfo, length_R_modeled - changepoint_min_distance)
      distance = as.integer(floor(changepoint_max_distance/spline_knot_distance))
      changepoint_min_distance = as.integer(ceiling(changepoint_min_distance/spline_knot_distance))
      if (distance < 2 || changepoint_min_distance < 1) {
        cli::cli_abort(paste(
          "The choosen `spline_knot_distance` is too long for the changepoint",
          "min and max distances. Please either choose a larger",
          "`changepoint_min_distance`/'changepoint_max_distance`",
          "or a smaller `spline_knot_distance`."
        ))
      }
      last_scp_knot <- sum(spline_knots$interior <= last_change)
      modeldata <- use_soft_changepoints(
        scp_length = length(spline_knots$interior),
        last_knot = last_scp_knot,
        distance = distance,
        min_distance = changepoint_min_distance,
        min_distance_tolerance = trend_change_tolerance,
        strictness_tol_k = strictness_tol_k,
        sharpness_boltzmann = sharpness_boltzmann,
        strictness_alpha = strictness_alpha,
        modeldata
      )
    },
    required = c(
      ".metainfo$length_R_modeled"
    ),
    modeldata = modeldata
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  # dummy data
  modeldata <- add_dummies_exponential_smoothing(modeldata)
  modeldata <- add_dummies_basis_splines2(modeldata)
  modeldata <- add_dummies_R_vari_selection(modeldata)
  modeldata <- add_dummies_R_vari(modeldata)
  modeldata <- add_dummies_smooth_derivative(modeldata)
  modeldata <- add_dummies_gp(modeldata)

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_changepoint_splines = c()
  )

  return(modeldata)
}

#'@title Estimate Rt with a smooth derivative
#'
#'@description This option estimates the effective reproduction number Rt over
#'  time as a smooth trend. It uses sparse smoothing splines placed on the
#'  first-order differences of Rt.
#'
#'@param spline_knot_distance Distance (in days) between spline knots for the
#'  penalized smoothing splines. Shorter distances increase flexibility of the
#'  Rt trajectory but also the daily Rt uncertainty.
#'@param trend_prior_shape Normal-Exponential-Gamma (NEG) prior (shape) for the
#'  Rt trend. The NEG prior is sparse, i.e. it has a strong peak at zero
#'  (transmission remains unchanged) and long tails (allows for occasional large
#'  positive or negative changes). Smaller shape parameters will lead to more
#'  sparseness (i.e. fewer but sharper changes in Rt). Note that when adjusting
#'  the shape, you will likely also have to adjust the scale.
#'@param trend_prior_scale Normal-Exponential-Gamma (NEG) prior (scale) for the
#'  strength of the Rt trend. See `sharp_changes_prior_shape` above for an
#'  explanation. Larger scales will lead to more Rt variability.
#'
#'@details The priors of this component have the following functional form:
#' - R_start (intercept):
#'  `Normal`
#'- trend_prior:
#'  `Normal-Exponential-Gamma`
#'
#'@inheritParams R_estimate_splines
#'
#'@inheritParams template_model_helpers
#'@inherit modeldata_init return
#'@export
#'@family {Rt models}
R_estimate_smooth_derivative <- function(
    R_start_prior_mu = 1,
    R_start_prior_sigma = 0.8,
    spline_knot_distance = 14,
    trend_prior_shape = 5,
    trend_prior_scale = 5e-3,
    link = "inv_softplus",
    R_max = 6,
    modeldata = modeldata_init()
) {

  modeldata <-  configure_R_model(
    name_approach = "smooth_derivative",
    model_id = 4,
    use_bs = TRUE,
    modeldata = modeldata
  )

  modeldata$R_intercept_prior <- set_prior(
    "R_intercept", "normal",
    mu = R_start_prior_mu, sigma = R_start_prior_sigma
  )
  modeldata$.init$R_intercept <-
    modeldata$R_intercept_prior$R_intercept_prior[1]

  modeldata$R_sd_change_prior <- set_prior("R_sd_change",
     "lomax",
     shape = trend_prior_shape,
     scale = trend_prior_scale
  )

  modeldata <- tbc(
    "R_smooth_derivative",
    {
      # Spline model
      spline_knots <- list(
        interior = rev(seq(
          modeldata$.metainfo$length_R_modeled - spline_knot_distance,
          1, by = -spline_knot_distance
        )),
        boundary = c(-1, modeldata$.metainfo$length_R_modeled)
      )
      modeldata$.metainfo$R_knots <- spline_knots
      modeldata <- use_basis_splines(
        spline_length = modeldata$.metainfo$length_R_modeled,
        knots = spline_knots,
        degree = 3,
        modeldata = modeldata
      )
      modeldata$.init$bs_coeff_noise_raw <- rep(1e-4, modeldata$bs_ncol - 2)
      modeldata$.init$bs_coeff_noise_lomax <- rep(1e-4, modeldata$bs_ncol - 2)
    },
    required = c(
      ".metainfo$length_R_modeled"
    ),
    modeldata = modeldata
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  # dummy data
  modeldata <- add_dummies_exponential_smoothing(modeldata)
  modeldata <- add_dummies_basis_splines2(modeldata)
  modeldata <- add_dummies_R_vari_selection(modeldata)
  modeldata <- add_dummies_R_vari(modeldata)
  modeldata <- add_dummies_soft_changepoints(modeldata)
  modeldata <- add_dummies_gp(modeldata)

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_smooth_derivative = c()
  )

  return(modeldata)
}

#' @title Estimate Rt using Gaussian processes
#'
#' @description This option estimates the effective reproduction number Rt over
#'   time using a Gaussian process (GP) model. There are two GPs: one for the
#'   long-term trend in Rt, and one for short-term deviations from this trend.
#'
#' @param R_intercept_prior_mu Prior (mean) for the intercept of Rt. Should be
#'   set to 1 unless you have a clear a priori expectation of the average Rt
#'   during the modeled time period.
#' @param R_intercept_prior_sigma Prior (standard deviation) for the intercept
#'   of Rt. By default, we fix `R_intercept` to 1 by setting
#'   `R_intercept_prior_sigma=0`. This is because deviations from Rt=1 are
#'   already captured by the long-term trend component (see below).
#' @param length_scale_prior_mu Prior (mean) on the length scale of the
#'   short-term Gaussian process (in days). This influences the smoothness of
#'   Rt. A higher length scale means that the Rt will change more slowly.
#'   Choosing a length scale that is too short can lead to overfitting of the Rt
#'   trajectory, while choosing a length scale that is too long can result in
#'   unrealistically smooth Rt estimates.
#' @param length_scale_prior_sigma Prior (standard deviation) on the length
#'   scale of the short-term Gaussian process. Set to zero to fix the length
#'   scale.
#' @param magnitude_prior_mu Prior (mean) on the magnitude of the short-term
#'   Gaussian process. Can be approximately interpreted as the marginal standard
#'   deviation of Rt. A higher magnitude allows more extreme Rt values.
#' @param magnitude_prior_sigma Prior (standard deviation) on the magnitude of
#'   the short-term Gaussian process. Set to zero to fix the magnitude.
#' @param long_length_scale_prior_mu Prior (mean) on the length scale of the
#'   long-term Gaussian process (in days). This should be quite long, at least
#'   several times the mean shedding delay of the pathogen, and significantly
#'   larger than the mean prior for the short-term GP (`length_scale_prior_mu`).
#' @param long_length_scale_prior_sigma Prior (standard deviation) on the length
#'   scale of the long-term Gaussian process. Set to zero to fix the length
#'   scale.
#' @param long_magnitude_prior_mu Prior (mean) on the magnitude of the long-term
#'   Gaussian process. Can be approximately interpreted as the marginal standard
#'   deviation of the long-term Rt trend. A higher magnitude allows more extreme
#'   Rt values.
#' @param long_magnitude_prior_sigma Prior (standard deviation) on the magnitude
#'   of the long-term Gaussian process. Set to zero to fix the magnitude.
#' @param matern_nu The smoothness parameter of the Matern kernel. The default
#'   is 3/2, other possible choices are 5/2 (more smooth) and 1/2 (less smooth).
#'   However, we recommend tuning smoothness primarily using the length scale
#'   priors.
#' @param boundary_factor The boundary factor used in the Gaussian process
#'   approximation. The default (`boundary_factor = 3`) is higher than the
#'   minimum recommendation from Riutort-Mayol et al. to ensure accurate
#'   real-time estimation. Lower values can lead to boundary effects close to
#'   the present and start of the time series and are thus not recommended. When
#'   modeling very long shedding delays, the boundary factor might need to be
#'   further increased. Note that this will also automatically increase the
#'   number of basis functions used and thereby slow down sampling.
#' @param n_basis_factor Factor used to automatically determine `m`, the number
#'   of basis functions in the Gaussian process approximation. Based on
#'   recommendations from Riutort-Mayol et al., we let `m = n_basis_factor * c /
#'   (l / S)`, where `c` is the `boundary_factor`, `l` is the length scale of
#'   the GP (we use the 5% quantile of `gp_length_prior` for a conservative
#'   result), and `S` is the maximum absolute value of the zero-centered input
#'   space, which is `(n-1)/2` where `n` is the number of Rt time steps modeled.
#'   This means that the number of basis function automatically adapts to the
#'   length scale and number of observations in the model. For the default
#'   settings (`boundary_factor = 3` and `n_basis_factor = 3.42`) and a lower
#'   length scale of 3 weeks, this corresponds to approx. 0.25x the number of Rt
#'   time points modeled. Increasing `n_basis_factor` will make the
#'   approximation more accurate by proportionally increasing the number of
#'   basis functions, but can slow down sampling.
#'
#' @details The estimated Rt trajectory is primarily influenced by the priors
#'   for the *length scale* and the *magnitude* of the short-term and long-term
#'   Gaussian processes. We recommend adjusting the *magnitude* when maximum Rt
#'   values seem to be too low or too high. In contrast, we recommend adjusting
#'   the *length scale* when the Rt trajectory seems to have too little
#'   resolution (try shorter length scales) or is overfitting on noise (try
#'   longer length scales). Note that, to ensure identifiability, the prior mean
#'   for the length scale of the long-term GP should always be significantly
#'   larger than the prior mean for the short-term GP.
#'
#' @details The Gaussian process is modeled using a Hilbert space approximation
#'   as described in Riutort-Mayol et al. (2023). This allows for fast inference
#'   without significant loss of accuracy. See the reference for more detail on
#'   choosing an adequate boundary factor and basis function factor. The
#'   `EpiSewer` implementation of the approximate Gaussian process is strongly
#'   inspired by the implementation in the
#'   [EpiNow2](https://github.com/epiforecasts/EpiNow2/) package.
#'
#' @details The priors of this component have the following functional form:
#' - R_intercept (intercept):
#'   `Normal`
#' - magnitude (magnitude):
#'   `Truncated Normal`
#' - length_scale (length scale):
#'   `Truncated Normal`
#' - long_magnitude (magnitude for long-term trend):
#'   `Truncated Normal`
#' - long_length_scale (length scale for long-term trend):
#'   `Truncated Normal`
#'
#' @references Riutort-Mayol, G., Bürkner, PC., Andersen, M.R. et al. Practical
#'   Hilbert space approximate Bayesian Gaussian processes for probabilistic
#'   programming. *Stat Comput* 33, 17 (2023).
#'   \url{https://doi.org/10.1007/s11222-022-10167-2}
#' @references Abbott S, Hellewell J, Thompson RN et al. Estimating the
#'   time-varying reproduction number of SARS-CoV-2 using national and
#'   subnational case counts. *Wellcome Open Res* 5, 112 (2020).
#'   \url{https://doi.org/10.12688/wellcomeopenres.16006.2}
#'
#' @inheritParams R_estimate_splines
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
R_estimate_gp <- function(
    R_intercept_prior_mu = 1,
    R_intercept_prior_sigma = 0,
    length_scale_prior_mu = 7*3,
    length_scale_prior_sigma = 7/2,
    magnitude_prior_mu = 0.2,
    magnitude_prior_sigma = 0.05,
    long_length_scale_prior_mu = 7*4*3,
    long_length_scale_prior_sigma = 7,
    long_magnitude_prior_mu = 0.4,
    long_magnitude_prior_sigma = 0.1,
    matern_nu = c(3/2, 5/2, 1/2),
    boundary_factor = 3,
    n_basis_factor = 3.42,
    link = "inv_softplus",
    R_max = 6,
    modeldata = modeldata_init()
) {

  modeldata <- configure_R_model(
    name_approach = "gp",
    model_id = 5,
    use_gp = TRUE,
    modeldata = modeldata
  )

  # intercept
  modeldata$R_intercept_prior <- set_prior(
    "R_intercept", "normal",
    mu = R_intercept_prior_mu, sigma = R_intercept_prior_sigma
  )
  modeldata$.init$R_intercept <- init_from_location_scale_prior(
    modeldata$R_intercept_prior
  )

  # magnitude
  modeldata$gp_sigma_prior <- set_prior("gp_sigma",
    mu = long_magnitude_prior_mu, sigma = long_magnitude_prior_sigma
  )
  modeldata$.init$gp_sigma <- init_from_location_scale_prior(
    modeldata$gp_sigma_prior
  )
  modeldata$gp2_sigma_prior <- set_prior("gp2_sigma",
    mu = magnitude_prior_mu, sigma = magnitude_prior_sigma
  )
  modeldata$.init$gp2_sigma <- init_from_location_scale_prior(
    modeldata$gp2_sigma_prior
  )

  # length scale
  modeldata$gp_length_prior <- set_prior("gp_length",
    mu = long_length_scale_prior_mu, sigma = long_length_scale_prior_sigma
  )
  modeldata$.init$gp_length <- init_from_location_scale_prior(
    modeldata$gp_length_prior
  )
  modeldata$gp_length_max <- long_length_scale_prior_mu + 3 * long_length_scale_prior_sigma
  modeldata$gp2_length_prior <- set_prior("gp2_length",
    mu = length_scale_prior_mu, sigma = length_scale_prior_sigma
  )
  modeldata$.init$gp2_length <- init_from_location_scale_prior(
    modeldata$gp2_length_prior
  )
  modeldata$gp2_length_max <- length_scale_prior_mu + 3 * length_scale_prior_sigma

  # further settings
  if (identical(matern_nu,c(3/2, 5/2, 1/2))) matern_nu <- 3/2
  if (matern_nu %in% c(3/2, 5/2, 1/2)) {
    modeldata$gp_matern_nu <- matern_nu
  } else {
    cli::cli_abort(
      "The Matern kernel parameter `nu` must be one of 3/2, 5/2, 1/2."
    )
  }
  modeldata$gp_c <- boundary_factor

  modeldata <- tbc(
    "R_gp",
    {
      modeldata$gp_n <- with(
        modeldata$.metainfo, length_R_modeled + forecast_horizon
        )
      S <- (modeldata$gp_n - 1) / 2.0 # maximum absolute value of input space

      l = max(1, long_length_scale_prior_mu - 2 * long_length_scale_prior_sigma) # length scale (5% quantile of prior)
      modeldata$gp_m <- ceiling(n_basis_factor * modeldata$gp_c / (l / S)) # number of basis functions
      modeldata$.init$gp_noise_raw <- rep(1e-4, modeldata$gp_m)

      l2 = max(1, length_scale_prior_mu - 2 * length_scale_prior_sigma) # length scale (5% quantile of prior)
      modeldata$gp2_m <- ceiling(n_basis_factor * modeldata$gp_c / (l2 / S)) # number of basis functions
      modeldata$.init$gp2_noise_raw <- rep(1e-4, modeldata$gp2_m)
    },
    required = c(
      ".metainfo$length_R_modeled",
      ".metainfo$forecast_horizon"
    ),
    modeldata = modeldata
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  # dummy data
  modeldata <- add_dummies_exponential_smoothing(modeldata)
  modeldata <- add_dummies_basis_splines(modeldata)
  modeldata <- add_dummies_basis_splines2(modeldata)
  modeldata <- add_dummies_R_vari_selection(modeldata)
  modeldata <- add_dummies_R_vari(modeldata)
  modeldata <- add_dummies_soft_changepoints(modeldata)
  modeldata <- add_dummies_smooth_derivative(modeldata)
  modeldata$R_sd_change_prior <- c(-1, -1)

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_gp = c()
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
add_R_variability <- function(length_R, h, length_seeding, partial_window,
                              partial_generation,
                              changepoint_dist, modeldata) {
  if (changepoint_dist == 0) {
    # no changepoints
    B <- matrix(rep(1,length_R+h), ncol = 1)
  } else if (length_R - partial_window - partial_generation - changepoint_dist <= length_seeding) {
    B <- matrix(rep(1,length_R+h), ncol = 1)
  } else {
    # we here suppress extrapolation warnings if h > changepoint_dist, as we fix
    # spline bases in the next step
    B <- suppressWarnings(splines::bs(
      1:(length_R+h),
      knots = rev(c(seq(
        length_R - partial_window - partial_generation,
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
  modeldata$R_vari_ncol <- ncol(B)
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
#'   [seeding_estimate_growth()] instead of [seeding_estimate_constant()].
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
    seeding_estimate_constant = c()
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
#'   detects (i.e. non-zero measurements). Before that date, the reproduction
#'   number will be retrospectively computed based on the seeded infections.
#'   Explicit modeling of the Rt time series will only begin after the seeding
#'   phase. This option avoids sampling problems when estimating Rt from very
#'   low infection numbers during a period with many non-detects. Note that
#'   estimated reproduction numbers are not necessarily meaningful during
#'   periods with very low infection numbers, as transmission dynamics may be
#'   dominated by chance events and importations.
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

#' Estimate seeding infections with a time-varying growth rate
#'
#' @description This option estimates an exponential growth of infections at the
#'   start of the modeled time period, with the growth rate varying over time.
#'
#' @param growth_change_prior_mu Prior (mean) on the daily standard deviation of
#'   the random walk for the epidemic growth rate.
#' @param growth_change_prior_sigma Prior (standard deviation) on the daily
#'   standard deviation of the random walk for the epidemic growth rate.
#'
#' @inheritParams seeding_estimate_rw
#'
#' @details The seeding phase has the length of the maximum generation time
#'   (during this time, the renewal model cannot be applied). It is here assumed
#'   that the expected number of new infections follows an exponential growth or
#'   decline process during this time period. The exponential growth rate of the
#'   seeding phase follows a random walk that ends with a growth rate
#'   representing the first estimated reproduction number at the start of the
#'   modeling phase. This means that the intercept of the seeding phase growth
#'   rate depends on the `R_start_prior` provided in the `R_estimate_*`
#'   component.
#'
#' @details If the lower and upper intervals of the prior for the initial number
#'   of infections, i.e. `intercept_prior_q5` or `intercept_prior_q95` are not
#'   specified by the user, `EpiSewer` will compute a rough median empirical
#'   estimate of the number of cases using the supplied wastewater measurements
#'   and shedding assumptions, and then infer the missing quantiles based on
#'   this. If none of the quantiles are provided, they are set to be roughly
#'   1/10 and 10 times the empirical median estimate. We note that this is a
#'   violation of Bayesian principles (data must not be used to inform priors) -
#'   but a neglectable one, since it only ensures that the seeding is modeled on
#'   the right order of magnitude and does not have relevant impacts on later Rt
#'   estimates.
#'
#' @details The priors of this component have the following functional form:
#' - initial number of infections (log scale): `Normal`
#' - standard deviation of the random walk on the growth rate: `Truncated normal`
#'
#'   The priors for these parameters are determined based on the user-supplied
#'   arguments, using appropriate transformations and the two-sigma-rule of
#'   thumb.
#'
#' @details Credits to Samuel Brand and the authors of the EpiAware toolkit for
#'   the idea to back-calculate the growth rate of the seeding phase from the
#'   initial reproduction number.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
seeding_estimate_growth <- function(
    intercept_prior_q5 = NULL,
    intercept_prior_q95 = NULL,
    growth_change_prior_mu = 0,
    growth_change_prior_sigma = 0.01,
    extend = TRUE,
    modeldata = modeldata_init()) {

  modeldata$seeding_model <- 2

  modeldata <- add_seeding_intercept_prior(
    intercept_prior_q5 = intercept_prior_q5,
    intercept_prior_q95 = intercept_prior_q95,
    calling_f_name = "seeding_estimate_growth",
    modeldata = modeldata
  )

  modeldata$iota_log_seed_sd_prior <- set_prior(
    "iota_log_seed_sd",
    "truncated normal",
    mu = growth_change_prior_mu,
    sigma = growth_change_prior_sigma
  )

  modeldata$.init$iota_log_seed_sd <- 1e-4
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
    seeding_estimate_growth = c()
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
      modeldata$.metainfo$infection_curve_crude$infections,
      ".metainfo$infection_curve_crude"
    )
    modeldata$.init$I_log <- tbe(
      log(modeldata$.metainfo$infection_curve_crude$infections),
      ".metainfo$infection_curve_crude"
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
