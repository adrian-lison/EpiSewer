
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


