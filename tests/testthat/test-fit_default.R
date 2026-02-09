
test_that("Default SARS-CoV-2 Zurich example works", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("Default Influenza A Zurich example works", {
  job <- EpiSewer(
    data = ww_data_influenza_Zurich,
    assumptions = ww_assumptions_influenza_Zurich,
    measurements = model_measurements(
      concentrations = concentrations_observe(replicate_col = "replicate_id"),
    ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("Basic forecast Influenza A Zurich example works", {
  job <- EpiSewer(
    data = ww_data_influenza_Zurich,
    assumptions = ww_assumptions_influenza_Zurich,
    measurements = model_measurements(
      concentrations = concentrations_observe(replicate_col = "replicate_id"),
    ),
    forecast = model_forecast(horizon = horizon_assume(7)),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_has_forecast(job))
}
)

test_that("Influenza A Zurich example works with complex model", {
  job <- EpiSewer(
    data = ww_data_influenza_Zurich,
    measurements = model_measurements(
      concentrations = concentrations_observe(replicate_col = "replicate_id"),
      noise = noise_estimate_dPCR(replicates = TRUE),
      LOD = LOD_estimate_dPCR()
    ),
    sampling = model_sampling(
      outliers = outliers_estimate(),
      sample_effects = sample_effects_estimate_weekday(
        effect_prior_mu = 0,
        effect_prior_sigma = 0.4
        )
    ),
    sewage = model_sewage(
      residence_dist = residence_dist_assume(residence_dist = c(0.8, 0.2)),
      flows = flows_observe()
    ),
    shedding = model_shedding(
      load_per_case = load_per_case_calibrate(min_cases = 15),
      load_variation = load_variation_estimate(),
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
      #shedding_dist = shedding_dist_assume(shedding_dist = ww_assumptions_influenza_Zurich$shedding_dist, shedding_reference = "symptom_onset"),
      incubation_dist = incubation_dist_assume(incubation_dist = c(1))
    ),
    infections = model_infections(
      generation_dist = generation_dist_assume(
        generation_dist = ww_assumptions_influenza_Zurich$generation_dist
        ),
      R = R_estimate_gp(),
      seeding = seeding_estimate_growth(),
      infection_noise = infection_noise_estimate(overdispersion = TRUE)
      ),
    forecast = model_forecast(
      horizon = horizon_assume(7),
      damping = damping_assume(damping = 0.9)
      ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_has_forecast(job))
}
)


