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
        modeldata$.metainfo$length_seeding <- length(generation_dist)
      },
      required_assumptions = "generation_dist",
      modeldata = modeldata
    )

    modeldata$.str$infections[["generation_dist"]] <- list(
      generation_dist_assume = c()
    )

    return(modeldata)
  }

#' Smooth Rt estimates via exponential smoothing
#'
#' @description This option estimates the effective reproduction number over
#'   time using exponential smoothing. It implements Holt's linear trend method
#'   with dampening through an innovations state space model with a level,
#'   trend, and dampening component.
#'
#' @param level_prior_mu Prior (mean) on the initial level of Rt.
#' @param level_prior_sigma Prior (standard deviation) on the initial level of
#'   Rt.
#' @param trend_prior_mu Prior (mean) on the initial trend of Rt.
#' @param trend_prior_sigma Prior (standard deviation) on the initial trend of
#'   Rt.
#' @param sd_prior_mu Prior (mean) on the standard deviation of the innovations.
#' @param sd_prior_sigma Prior (standard deviation) on the standard deviation of
#'   the innovations.
#' @param smooth_fixed Fixed value for the smoothing parameter
#'   (alpha). If NULL (default), the smoothing parameter is estimated, otherwise
#'   it will be fixed to the given value.
#' @param smooth_prior_shapes Prior (Beta shape 1 and 2) on the smoothing
#'   parameter.
#' @param trend_smooth_fixed Fixed value for the trend smoothing parameter
#'   (beta). If NULL (default), the trend smoothing parameter is estimated,
#'   otherwise it will be fixed to the given value.
#' @param trend_smooth_prior_shapes Prior (Beta shape 1 and 2) on the trend
#'   smoothing parameter.
#' @param dampen_fixed Fixed value for the dampening parameter (phi). If NULL
#'   (default), the dampening parameter is estimated, otherwise it will be fixed
#'   to the given value.
#' @param dampen_prior_shapes Prior (Beta shape 1 and 2) on the dampening
#'   parameter.
#' @param differenced If `FALSE` (default), exponential smoothing is applied to
#'   the absolute Rt time series. If `TRUE`, it is instead applied to the
#'   differenced time series. This makes the level become the trend, and the
#'   trend become the curvature.
#' @param noncentered If `TRUE` (default), a non-centered parameterization is
#'   used to model the innovations in the state space process (for better
#'   sampling efficiency).
#'
#' @inheritParams R_estimate_splines
#'
#' @details The innovations state space model consists of three components: a
#'   level, a trend, and a dampening component.
#' - The level is smoothed based on the levels from earlier time steps,
#'   with exponentially decaying weights, as controlled by a smoothing parameter
#'   (often called alpha). Note that *smaller* values of `alpha` indicate
#'   *stronger* smoothing. In particular, `alpha = 1` means that only the last
#'   level is used.
#' - The trend is smoothed based on the trends from earlier time steps,
#'   with exponentially decaying weights, as controlled by a trend smoothing
#'   parameter (often called beta). Note that *smaller* values of `beta`
#'   indicate *stronger* smoothing. In particular, `beta = 1` means that only
#'   the last trend is used.
#' - The dampening determines how long a previous trend continues into the
#'   future before it levels of to a stationary time series. The strength of
#'   dampening is controlled by a dampening parameter (often called phi).
#'   Note that *smaller* values of `phi` indicate *stronger* dampening.
#'   In particular, `phi = 1` means no dampening. Values below `phi = 0.8`
#'   are seldom in practice as the dampening becomes very strong.
#'
#' @details Often, `alpha`, `beta`, and `phi` are jointly unidentifiable. It may
#'   therefore be necessary to fix at least one of the parameters (typically
#'   `phi`) or supply strong priors.
#'
#' @details Note that the smoothness of retrospective Rt estimates is often more
#'   influenced by the prior on the standard deviation of innovations than the
#'   smoothing and trend smoothing parameters. The smoothing parameters mostly
#'   have an influence on estimates close to the present / date of estimation,
#'   when limited data signal is available. Here, the standard deviation of the
#'   innovations influences how uncertain Rt estimates are close to the present.
#'
#' @details The priors of this component have the following functional form:
#' - initial level of Rt: `Normal`
#' - initial trend of Rt: `Normal`
#' - standard deviation of innovations: `Truncated normal`
#' - smoothing parameter: `Beta`
#' - trend smoothing parameter: `Beta`
#' - dampening parameter: `Beta`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
R_estimate_ets <- function(
    level_prior_mu = 1,
    level_prior_sigma = 0.8,
    trend_prior_mu = 0,
    trend_prior_sigma = 0.1,
    sd_prior_mu = 0,
    sd_prior_sigma = 0.1,
    link = "inv_softplus",
    R_max = 6,
    smooth_fixed = NULL,
    smooth_prior_shapes = c(50, 50),
    trend_smooth_fixed = NULL,
    trend_smooth_prior_shapes = c(50, 50),
    dampen_fixed = 0.9,
    dampen_prior_shapes = c(50, 5),
    differenced = FALSE,
    noncentered = TRUE,
    modeldata = modeldata_init()) {
  modeldata$.metainfo$R_estimate_approach <- "ets"

  modeldata$R_level_start_prior <- set_prior(
    "R_level_start", "normal",
    mu = level_prior_mu, sigma = level_prior_sigma
  )

  modeldata$R_trend_start_prior <- set_prior(
    "R_trend_start", "normal",
    mu = trend_prior_mu, sigma = trend_prior_sigma
  )

  modeldata$R_sd_prior <- set_prior(
    "R_sd", "truncated normal",
    mu = sd_prior_mu, sigma = sd_prior_sigma
  )

  modeldata$.init$R_level_start <-
    modeldata$R_level_start_prior$R_level_start_prior[1]
  modeldata$.init$R_trend_start <- 1e-4
  modeldata$.init$R_sd <- max(modeldata$R_sd_prior$R_sd_prior[1], 0.1)
  modeldata$.init$R_noise <- tbe(
    rep(0, modeldata$.metainfo$length_R - 1),
    ".metainfo$length_R"
  )

  if (is.null(smooth_fixed)) {
    smooth_fixed <- -1
  }
  modeldata$ets_alpha_fixed <- smooth_fixed
  if (smooth_fixed >= 0) {
    modeldata$ets_alpha_prior <- numeric(0)
    modeldata$.init$ets_alpha <- numeric(0)
  } else {
    modeldata$ets_alpha_prior <- set_prior("ets_alpha", "beta",
      alpha = smooth_prior_shapes[1], beta = smooth_prior_shapes[2]
    )
    modeldata$.init$ets_alpha <- as.array(0.5)
  }

  if (is.null(trend_smooth_fixed)) {
    trend_smooth_fixed <- -1
  }
  modeldata$ets_beta_fixed <- trend_smooth_fixed
  if (trend_smooth_fixed >= 0) {
    modeldata$ets_beta_prior <- numeric(0)
    modeldata$.init$ets_beta <- numeric(0)
  } else {
    modeldata$ets_beta_prior <- set_prior("ets_beta", "beta",
      alpha = trend_smooth_prior_shapes[1], beta = trend_smooth_prior_shapes[2]
    )
    modeldata$.init$ets_beta <- as.array(0.5)
  }

  if (is.null(dampen_fixed)) {
    dampen_fixed <- -1
  }
  modeldata$ets_phi_fixed <- dampen_fixed
  if (dampen_fixed >= 0) {
    modeldata$ets_phi_prior <- numeric(0)
    modeldata$.init$ets_phi <- numeric(0)
  } else {
    modeldata$ets_phi_prior <- set_prior("ets_phi", "beta",
      alpha = dampen_prior_shapes[1], beta = dampen_prior_shapes[2]
    )
    modeldata$.init$ets_phi <- as.array(0.9)
  }

  modeldata <- add_link_function(link, R_max, modeldata)

  modeldata$ets_diff <- differenced
  modeldata$ets_noncentered <- noncentered

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_ets = c()
  )

  return(modeldata)
}

#' Smooth Rt estimates via a random walk
#'
#' @description This option estimates the effective reproduction number over
#'   time using a random walk.
#'
#' @param intercept_prior_mu Prior (mean) on the intercept of the random walk.
#' @param intercept_prior_sigma Prior (standard deviation) on the intercept of
#'   the random walk.
#' @param sd_prior_mu Prior (mean) on the standard deviation of the random walk.
#' @param sd_prior_sigma Prior (standard deviation) on the standard deviation of
#'   the random walk.
#' @param differenced If `FALSE` (default), the random walk is applied to the
#'   absolute Rt time series. If `TRUE`, it is instead applied to the
#'   differenced time series, i.e. now the trend is modeled as a random walk.
#' @param noncentered If `TRUE` (default), a non-centered parameterization is
#'   used to model the innovations of the random walk (for better sampling
#'   efficiency).
#' @inheritParams R_estimate_ets
#'
#' @details The smoothness of Rt estimates is influenced by the prior on the
#'   standard deviation of the random walk. It also influences the uncertainty
#'   of Rt estimates towards the present / date of estimation, when limited
#'   data signal is available. The prior on the intercept of the random walk
#'   should reflect your expectation of Rt at the beginning of the time series.
#'   If estimating from the start of an epidemic, you might want to use a prior
#'   with mean > 1 for the intercept.
#'
#' @details The priors of this component have the following functional form:
#' - intercept of the random walk: `Normal`
#' - standard deviation of the random walk: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
R_estimate_rw <- function(
    intercept_prior_mu = 1,
    intercept_prior_sigma = 0.8,
    sd_prior_mu = 0,
    sd_prior_sigma = 0.1,
    link = "inv_softplus",
    R_max = 6,
    differenced = FALSE,
    noncentered = TRUE,
    modeldata = modeldata_init()) {
  modeldata <- R_estimate_ets(
    level_prior_mu = intercept_prior_mu,
    level_prior_sigma = intercept_prior_sigma,
    sd_prior_mu = sd_prior_mu,
    sd_prior_sigma = sd_prior_sigma,
    link = link,
    R_max = R_max,
    smooth_fixed = 1,
    trend_smooth_fixed = 0,
    dampen_fixed = 0,
    differenced = differenced,
    noncentered = noncentered,
    modeldata = modeldata
  )

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_rw = c()
  )

  return(modeldata)
}

#' Smooth Rt estimates via smoothing splines
#'
#' @description This option estimates the effective reproduction number using
#'   penalized B-splines.
#'
#' @param knot_distance Distance between spline breakpoints (knots) in days
#'   (default is 7, i.e. one knot each week). Placing knots further apart
#'   increases the smoothness of Rt estimates and can speed up model fitting.
#'   The Rt time series remains surprisingly flexible even at larger knot
#'   distances, but placing knots too far apart can lead to inaccurate
#'   estimates.
#' @param spline_degree Degree of the spline polynomials (default is 3 for cubic
#'   splines).
#' @param coef_intercept_prior_mu Prior (mean) on the intercept of the random
#'   walk over spline coefficients.
#' @param coef_intercept_prior_sigma Prior (standard deviation) on the intercept
#'   of the random walk over spline coefficients.
#' @param coef_sd_prior_mu Prior (mean) on the daily standard deviation of the
#'   random walk over spline coefficients.
#' @param coef_sd_prior_sigma Prior (standard deviation) on the daily standard
#'   deviation of the random walk over spline coefficients.
#' @param link Link function. Currently supported are `inv_softplus` (default)
#'   and `scaled_logit`. Both of these links are configured to behave
#'   approximately like the identity function around R=1, but become
#'   increasingly non-linear below (and in the case of `scaled_logit` also
#'   above) R=1.
#' @param R_max If `link=scaled_logit` is used, a maximum reproduction number
#'   must be assumed. This should be higher than any realistic R value for the
#'   modeled pathogen. Default is 6.
#'
#' @details `EpiSewer` uses a random walk on the B-spline coefficients for
#'   regularization. This allows to use small knot distances without obtaining
#'   extremely wiggly Rt estimates.
#'   - The prior on the  random walk intercept should reflect
#'   your expectation of Rt at the beginning of the time series. If estimating
#'   from the start of an epidemic, you might want to use a prior with mean
#'   larger than 1 for the intercept.
#'   - The prior on the daily standard deviation should be interpreted in terms of
#'   daily additive changes (this is at least true around Rt=1, and becomes less
#'   true as Rt approaches 0 or its upper bound as defined by the `link`
#'   function). For example, a standard deviation of 0.5 roughly allows the
#'   spline coefficients to change by ±1 (using the 2 sigma rule) each day. The
#'   daily standard deviation is multiplied by `sqrt(knot_distance)` to get the
#'   actual standard deviation of the coefficients. This way, the prior is
#'   independent of the chosen `knot_distance`. For example, if `knot_distance`
#'   is 7 days, and a daily standard deviation of 0.5 is used, the coefficients
#'   of two adjacent splines can differ by up to ±2.6. Note however that this
#'   difference does not directly translate into a change of Rt by ±2.6, as Rt
#'   is always the weighted sum of several basis functions at any given point.
#'   It may therefore change much more gradually, depending on the distances
#'   between knots.
#'
#' @details The smoothness of the Rt estimates is influenced both by the knot
#'   distance and by the daily standard deviation of the random walk on
#'   coefficients. The latter also influences the uncertainty of Rt estimates
#'   towards the present / date of estimation, when limited data signal is
#'   available.
#'
#' @details The priors of this component have the following functional form:
#' - intercept of the log random walk over spline coefficients:
#'   `Normal`
#' - standard deviation of the random walk over spline coefficients:
#'   `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {Rt models}
R_estimate_splines <- function(
    knot_distance = 7,
    spline_degree = 3,
    coef_intercept_prior_mu = 1,
    coef_intercept_prior_sigma = 0.8,
    coef_sd_prior_mu = 0,
    coef_sd_prior_sigma = 0.5,
    link = "inv_softplus",
    R_max = 6,
    modeldata = modeldata_init()) {
  modeldata$.metainfo$R_estimate_approach <- "splines"

  modeldata <- tbc(
    "spline_definition",
    {
      knots <- place_knots(
        ts_length = modeldata$.metainfo$length_R,
        knot_distance = knot_distance,
        partial_window = modeldata$.metainfo$partial_window
      )
      B <-
        splines::bs(
          1:modeldata$.metainfo$length_R,
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

      modeldata$.init$bs_coeff_ar_start <- 0
      modeldata$.init$bs_coeff_ar_sd <- 0.1
      modeldata$.init$bs_coeff_noise_raw <- rep(0, modeldata$bs_n_basis - 1)
    },
    required = c(".metainfo$length_R", ".metainfo$partial_window"),
    modeldata = modeldata
  )

  modeldata$bs_coeff_ar_start_prior <- set_prior("bs_coeff_ar_start",
    "normal",
    mu = coef_intercept_prior_mu,
    sigma = coef_intercept_prior_sigma
  )
  modeldata$bs_coeff_ar_sd_prior <- set_prior("bs_coeff_ar_sd",
    "truncated normal",
    mu = coef_sd_prior_mu,
    sigma = coef_sd_prior_sigma
  )

  modeldata <- add_link_function(link, R_max, modeldata)

  modeldata$.str$infections[["R"]] <- list(
    R_estimate_splines = c()
  )

  return(modeldata)
}

#' Estimate seeding infections
#'
#' @description This option estimates initial infections at the start of the
#'   modeled time period, when the renewal model cannot be applied yet. It uses
#'   a random walk to model these seeding infections.
#'
#' @param intercept_prior_mu Prior (mean) on the intercept of the seeding random
#'   walk. If NULL (default), this is set to a crude empirical estimate of the
#'   number of cases (see details).
#' @param intercept_prior_sigma Prior (standard deviation) on the intercept of
#'   the seeding random walk.
#' @param sd_prior_mu Prior (mean) on the standard deviation of the seeding
#'   random walk.
#' @param sd_prior_sigma Prior (standard deviation) on the standard deviation of
#'   the seeding random walk.
#'
#' @details The seeding phase has the length of the maximum generation time
#'   (during this time, the renewal model cannot be applied). Traditionally,
#'   seeding refers to the first few (potentially imported) infections of an
#'   epidemic, but depending on what time period the model is fitted to, this
#'   may also cover a different phase with higher incidence levels.
#'
#' @details If `intercept_prior_mu` is not specified by the user, `EpiSewer`
#'   will set it to a rough estimate of the number of cases using the supplied
#'   wastewater measurements and shedding assumptions. We note that this is a
#'   violation of Bayesian principles (data must not be used to inform priors) -
#'   but a neglectable one, since it only ensures that the seeding is modeled on
#'   the right order of magnitude and does not have relevant impacts on later Rt
#'   estimates.
#'
#' @details The priors of this component have the following functional form:
#' - intercept of the seeding random walk: `Normal`
#' - standard deviation of the seeding random walk: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
seeding_estimate_rw <- function(
    intercept_prior_mu = NULL,
    intercept_prior_sigma = 1,
    sd_prior_mu = 0.05,
    sd_prior_sigma = 0.025,
    modeldata = modeldata_init()) {
  if (!is.null(intercept_prior_mu)) {
    modeldata$iota_log_ar_start_prior <- set_prior(
      "iota_log_ar_start",
      "normal",
      mu = intercept_prior_mu,
      sigma = intercept_prior_sigma
    )
  } else {
    modeldata$iota_log_ar_start_prior <- tbe(
      set_prior(
        "iota_log_ar_start",
        "normal (mu based on crude empirical estimate of cases)",
        mu = log(modeldata$.metainfo$initial_cases_crude),
        sigma = intercept_prior_sigma
      ),
      ".metainfo$initial_cases_crude"
    )
  }

  modeldata$iota_log_ar_sd_prior <- set_prior(
    "iota_log_ar_sd",
    "truncated normal",
    mu = sd_prior_mu,
    sigma = sd_prior_sigma
  )

  modeldata$.init$iota_log_ar_start <- tbe(
    log(modeldata$.metainfo$initial_cases_crude),
    ".metainfo$initial_cases_crude"
  )
  modeldata$.init$iota_log_ar_sd <- 1
  modeldata$.init$iota_log_ar_noise <- tbe(
    rep(0, modeldata$.metainfo$length_seeding - 1),
    ".metainfo$length_seeding"
  )

  modeldata$.str$infections[["seeding"]] <- list(
    seeding_estimate_rw = c()
  )

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
#' @param overdispersion If `FALSE` (default) new infections are modeled as
#'   Poisson distributed. If `TRUE`, new infections are modeled as Negative
#'   Binomial distributed.
#' @param overdispersion_prior_mu Prior (mean) on the overdispersion parameter
#'   of the Negative Binomial.
#' @param overdispersion_prior_sigma Prior (standard deviation) on the
#'   overdispersion parameter of the Negative Binomial.
#'
#' @details The priors of this component have the following functional form:
#' - overdispersion parameter of the Negative Binomial: `Truncated normal`
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {infection noise models}
infection_noise_estimate <-
  function(overdispersion = FALSE,
           overdispersion_prior_mu = 0,
           overdispersion_prior_sigma = 1,
           modeldata = modeldata_init()) {
    modeldata$I_overdispersion <- overdispersion
    modeldata$I_xi_prior <- set_prior("I_xi", "normal",
      mu = overdispersion_prior_mu, sigma = overdispersion_prior_sigma
    )

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

    if (modeldata$I_overdispersion) {
      modeldata$.init$I_xi <- as.array(0.05)
    } else {
      modeldata$I_xi_prior <- numeric(0)
      modeldata$.init$I_xi <- numeric(0)
    }

    modeldata$.str$infections[["infection_noise"]] <- list(
      infection_noise_estimate = c()
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
      rlang::abort("Please chose a maximum R value greater or equal to 3.")
    }
    modeldata$.metainfo$R_max <- R_max
    a_k <- logistic_find_a_k(R_max)
    # first element: 1 = scaled_logit
    # hyperparameters are chosen such that R=1 at x=1, and is roughly
    # linear with slope 1 around R=1
    modeldata$R_link <- c(1, R_max, a_k$a, a_k$k)
  } else {
    rlang::abort(c(
      paste0("'", link, "' is not a valid link function. ",
             "Please use one of:"), "inv_softplus", "scaled_logit"))
  }
  return(modeldata)
}
