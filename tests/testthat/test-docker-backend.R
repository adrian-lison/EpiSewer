
test_that("Running from R with docker backend works", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    fit_opts = set_fit_opts(
      sampler = sampler_stan_mcmc(
        iter_warmup = 5, iter_sampling = 5, chains = 2
      ),
      model = model_stan_opts(use_docker = TRUE)
    ),
    run_fit = FALSE
  )
  job <- run(job, run_silent = TRUE)
  expect_true(check_result_valid(job))
}
)
