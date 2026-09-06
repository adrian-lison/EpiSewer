# Battery of model configurations for modeldata characterization tests.
#
# Each entry is a zero-argument function returning a list of arguments for
# EpiSewer(). The same definitions drive both the one-off fixture recording
# script (data-raw/record_modeldata_fixtures.R) and the comparison tests
# (test-modeldata_characterization.R), so that recorded and rebuilt jobs
# always come from identical specifications.
#
# Together the configs cover every model component helper at least once.

modeldata_config_battery <- function() {
  list(

    # 1: all defaults (concentrations_observe, noise_estimate, LOD_none,
    # outliers_estimate, sample_effects_none, flows_observe, residence default,
    # shedding/incubation assume via symptom_onset, load_per_case_calibrate
    # via cases data, R_estimate_gp, seeding_estimate_rw,
    # infection_noise_estimate, horizon_none, damping default)
    default_sars = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich
      )
    },

    # 2: dPCR noise with replicates
    dpcr_replicates = function() {
      list(
        data = ww_data_influenza_Zurich,
        assumptions = ww_assumptions_influenza_Zurich,
        measurements = model_measurements(
          concentrations = concentrations_observe(replicate_col = "replicate_id"),
          noise = noise_estimate_dPCR(replicates = TRUE),
          LOD = LOD_estimate_dPCR()
        )
      )
    },

    # 3: dPCR noise with estimated assay params, total partitions observed
    dpcr_params_partitions_observed = function() {
      list(
        data = ww_data_influenza_Zurich,
        assumptions = ww_assumptions_influenza_Zurich,
        measurements = model_measurements(
          concentrations = concentrations_observe(
            replicate_col = "replicate_id",
            total_partitions_col = "total_droplets",
            n_averaged = 2
          ),
          noise = noise_estimate_dPCR_params(total_partitions_observe = TRUE),
          LOD = LOD_estimate_dPCR()
        )
      )
    },

    # 4: dPCR noise with estimated assay params, partition number estimated
    # (max_partitions / partition_loss prior branch)
    dpcr_params_estimated = function() {
      list(
        data = ww_data_influenza_Zurich,
        assumptions = ww_assumptions_influenza_Zurich,
        measurements = model_measurements(
          concentrations = concentrations_observe(replicate_col = "replicate_id"),
          noise = noise_estimate_dPCR_params(),
          LOD = LOD_estimate_dPCR()
        )
      )
    },

    # 5: binomial likelihood on positive partition counts
    partitions_likelihood = function() {
      list(
        data = ww_data_influenza_Zurich,
        assumptions = ww_assumptions_influenza_Zurich,
        measurements = model_measurements(
          concentrations = concentrations_observe_partitions(
            positive_partitions_col = "positive_droplets",
            total_partitions_col = "total_droplets",
            replicate_col = "replicate_id"
          ),
          noise = noise_estimate_dPCR(total_partitions_observe = TRUE),
          LOD = LOD_estimate_dPCR()
        )
      )
    },

    # 6: assumed limit of detection
    lod_assume = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        measurements = model_measurements(
          LOD = LOD_assume(limit = 2.5, prob = 0.95)
        )
      )
    },

    # 7: constant variance noise, log-normal observation distribution
    constant_var = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        measurements = model_measurements(
          noise = noise_estimate_constant_var(
            distribution = "log-normal", warn = FALSE
          )
        )
      )
    },

    # 8: ETS Rt model with forecast, no damping
    ets_forecast = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        infections = model_infections(R = R_estimate_ets()),
        forecast = model_forecast(
          horizon = horizon_assume(7),
          damping = damping_none()
        )
      )
    },

    # 9: random walk Rt model
    rw = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        infections = model_infections(R = R_estimate_rw())
      )
    },

    # 10: spline Rt model
    splines = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        infections = model_infections(R = R_estimate_splines())
      )
    },

    # 11: changepoint spline Rt model with damped forecast
    changepoint_splines = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        infections = model_infections(R = R_estimate_changepoint_splines()),
        forecast = model_forecast(
          horizon = horizon_assume(14),
          damping = damping_assume(damping = 0.9)
        )
      )
    },

    # 12: smooth derivative Rt model
    smooth_derivative = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        infections = model_infections(R = R_estimate_smooth_derivative())
      )
    },

    # 13: piecewise Rt model
    piecewise = function() {
      list(
        data = ww_data_SARS_CoV_2_Zurich,
        assumptions = ww_assumptions_SARS_CoV_2_Zurich,
        infections = model_infections(R = R_estimate_piecewise())
      )
    },

    # 14: alternative components: assumed flows, weekday sample effects, no
    # outliers, constant seeding, no infection noise, load variation model,
    # residence distribution, shedding referenced by infection
    alt_components = function() {
      assumptions <- ww_assumptions_SARS_CoV_2_Zurich
      assumptions$shedding_reference <- "infection"
      list(
        data = sewer_data(
          measurements = ww_data_SARS_CoV_2_Zurich$measurements
        ),
        assumptions = assumptions,
        sampling = model_sampling(
          outliers = outliers_none(),
          sample_effects = sample_effects_estimate_weekday()
        ),
        sewage = model_sewage(
          flows = flows_assume(flow_constant = 1e8),
          residence_dist = residence_dist_assume(residence_dist = c(0.8, 0.2))
        ),
        shedding = model_shedding(
          shedding_dist = shedding_dist_assume(
            shedding_dist = ww_assumptions_SARS_CoV_2_Zurich$shedding_dist,
            shedding_reference = "infection"
          ),
          load_variation = load_variation_estimate()
        ),
        infections = model_infections(
          seeding = seeding_estimate_constant(),
          infection_noise = infection_noise_none()
        )
      )
    },

    # 15: estimated shedding distribution (mixture), growth seeding, assumed
    # load per case, design-matrix sample effects, forecast
    estimate_dists_growth = function() {
      assumptions <- ww_assumptions_influenza_Zurich
      assumptions$shedding_reference <- "symptom_onset"
      assumptions$shedding_dist <- NULL
      measurements <- ww_data_influenza_Zurich$measurements
      T_days <- as.integer(
        max(measurements$date, na.rm = TRUE) -
          min(measurements$date, na.rm = TRUE)
      ) + 1L
      h <- 7L
      day_index <- seq_len(T_days + h)
      design_matrix <- cbind(
        sinus = sin(2 * pi * day_index / 7),
        cosinus = cos(2 * pi * day_index / 7)
      )
      list(
        data = ww_data_influenza_Zurich,
        assumptions = assumptions,
        measurements = model_measurements(
          concentrations = concentrations_observe(replicate_col = "replicate_id")
        ),
        sampling = model_sampling(
          sample_effects = sample_effects_estimate_matrix(
            design_matrix = design_matrix
          )
        ),
        shedding = model_shedding(
          shedding_dist = shedding_dist_estimate(
            shedding_dist_mean_prior_mean = c(0.1, 0.12),
            shedding_dist_mean_prior_sd = c(0.002, 0.001),
            shedding_dist_cv_prior_mean = c(0.5, 0.6),
            shedding_dist_cv_prior_sd = c(0.02, 0.01),
            shedding_dist_type = "gamma",
            shedding_reference = "symptom_onset",
            prior_weights = c(0.7, 0.3),
            weight_alpha = 1
          ),
          incubation_dist = incubation_dist_assume(incubation_dist = c(1)),
          load_per_case = load_per_case_assume(load_per_case = 1e11)
        ),
        infections = model_infections(
          seeding = seeding_estimate_growth()
        ),
        forecast = model_forecast(horizon = horizon_assume(h))
      )
    },

    # 16: calibration to case data, data passed directly to helpers,
    # leading non-detects triggering the seeding extension
    calibrate_cases_direct = function() {
      measurements <- data.table::copy(ww_data_influenza_Zurich$measurements)
      zero_until <- min(measurements$date, na.rm = TRUE) + 13
      measurements[
        date <= zero_until & !is.na(concentration),
        c("concentration", "positive_droplets") := list(0, 0)
      ]
      flows <- ww_data_influenza_Zurich$flows
      cases <- data.frame(
        date = sort(unique(flows$date)),
        cases = 25 + 10 * sin(seq_along(unique(flows$date)) / 20)
      )
      list(
        assumptions = ww_assumptions_influenza_Zurich,
        measurements = model_measurements(
          concentrations = concentrations_observe(
            measurements = measurements,
            replicate_col = "replicate_id"
          )
        ),
        sewage = model_sewage(
          flows = flows_observe(flows = flows)
        ),
        shedding = model_shedding(
          load_per_case = load_per_case_calibrate(cases = cases)
        )
      )
    }
  )
}

# Fields of job$metainfo that are compared against fixtures (stable external
# contract; excludes crude-curve tables and R model internals).
modeldata_fixture_meta_fields <- function() {
  c(
    "T_start_date", "T_end_date", "forecast_horizon", "measured_dates",
    "total_delay_dist", "measurements_cols", "composite_window",
    "load_per_case", "length_seeding", "length_R", "length_R_modeled",
    "length_I", "date_triple_detect"
  )
}

# Build the (unfitted) EpiSewer job for one battery config.
# Names are sorted with a locale-independent radix sort (testthat runs under
# LC_COLLATE=C, the recording script under the user locale).
modeldata_fixture_job <- function(config_args) {
  sort_c <- function(x) sort(x, method = "radix")
  job <- do.call(EpiSewer, c(config_args, list(run_fit = FALSE)))$job
  meta <- job$metainfo[
    intersect(modeldata_fixture_meta_fields(), names(job$metainfo))
  ]
  list(
    data = job$data[sort_c(names(job$data))],
    init = job$init[sort_c(names(job$init))],
    metainfo = meta[sort_c(names(meta))]
  )
}

modeldata_fixture_path <- function(name) {
  testthat::test_path("fixtures", "modeldata", paste0(name, ".rds"))
}
