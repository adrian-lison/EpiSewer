dist1 <- get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397)
dist2 <- get_discrete_gamma(gamma_shape = 0.429639, gamma_scale = 7.241397)
dist3 <- get_discrete_gamma(gamma_shape = 3.429639, gamma_scale = 3.241397)

test_that("shedding_dist_estimate_mixture works", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    shedding = model_shedding(shedding_dist = shedding_dist_estimate_mixture(
      shedding_dist_list = list(dist1, dist2, dist3),
      prior_weights = c(0.2, 0.5, 0.3),
    )),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("shedding_dist_estimate_mixture works different weights and alpha", {
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_assumptions_SARS_CoV_2_Zurich,
    shedding = model_shedding(shedding_dist = shedding_dist_estimate_mixture(
      shedding_dist_list = list(dist1, dist2, dist3),
      prior_weights = c(5, 2, 4), weight_alpha = 0.6
    )),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)

test_that("shedding_dist_estimate_mixture works with reference to infection ", {
  ww_ass <- ww_assumptions_SARS_CoV_2_Zurich
  ww_ass$shedding_reference <- "infection"
  job <- EpiSewer(
    data = ww_data_SARS_CoV_2_Zurich,
    assumptions = ww_ass,
    shedding = model_shedding(shedding_dist = shedding_dist_estimate_mixture(
      shedding_dist_list = list(dist1, dist2, dist3),
      prior_weights = c(0.2, 0.5, 0.3),
    )),
    run_fit = FALSE
  )
  job <- suppressWarnings(test_run(job))
  expect_true(check_result_valid(job))
}
)
