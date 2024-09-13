test_that("EpiSewer fails are detected", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    run_fit = FALSE
  )
  job$job$data$L = -1
  job <- suppressWarnings(test_run(job))
  expect_error(check_result_valid(job))
}
)
