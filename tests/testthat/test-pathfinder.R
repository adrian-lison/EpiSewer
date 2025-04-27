
test_that("Default model can be fitted with Pathfinder", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    run_fit = FALSE,
    fit_opts = set_fit_opts(
      sampler = sampler_stan_pathfinder()
    )
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("Default model can be fitted with Pathfinder initialization", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    run_fit = FALSE,
    fit_opts = set_fit_opts(
      sampler = sampler_stan_mcmc(init_pathfinder = TRUE)
    )
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)
