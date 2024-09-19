test_that("Spline smoothing with 1 degrees of freedom works", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    infections = model_infections(
      R = R_estimate_splines(robust = TRUE, robust_df = 1)
    ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("Spline smoothing with 4 degrees of freedom works", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    infections = model_infections(
      R = R_estimate_splines(robust = TRUE, robust_df = 4)
    ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("Spline smoothing with infinite degrees of freedom works", {
  # infinite df equals to Gaussian errors instead of Student-t
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    infections = model_infections(
      R = R_estimate_splines(robust = FALSE)
    ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("ETS smoothing with 1 degrees of freedom works", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    infections = model_infections(
      R = R_estimate_ets(robust = TRUE, robust_df = 1)
    ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("ETS smoothing with 4 degrees of freedom works", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    infections = model_infections(
      R = R_estimate_ets(robust = TRUE, robust_df = 4)
    ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("ETS smoothing with infinite degrees of freedom works", {
  # infinite df equals to Gaussian errors instead of Student-t
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    infections = model_infections(
      R = R_estimate_ets(robust = FALSE)
    ),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)


