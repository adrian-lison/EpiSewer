md <- concentrations_observe(
  measurements = ww_data_SARS_CoV_2_Zurich$measurements[date <= "2022-04-07",]
  )
# set forecasting horizon
md <- horizon_assume(horizon = 7, modeldata = md)

test_that("flows_observe() throws error when flow is missing on sampled dates", {
  # take 5 dates from the sampled dates
  missing_dates <- ww_data_SARS_CoV_2_Zurich$measurements$date[c(1, 5, 7, 8, 12)]
  flow_missing <- ww_data_SARS_CoV_2_Zurich$flows[!(date %in% missing_dates),]
  expect_error(
    flows_observe(flows = flow_missing, modeldata = md),
    "Missing flow values for the following sampled dates: 2022-01-03,"
    )
}
)

test_that("flows_observe() imputes flow values missing for unsampled dates", {
  # take 5 unsampled dates
  missing_dates <- 1 + ww_data_SARS_CoV_2_Zurich$measurements$date[c(1, 5, 7, 8, 12)]
  flow_missing <- ww_data_SARS_CoV_2_Zurich$flows[!(date %in% missing_dates),]
  expect_message(
    md_flow <- flows_observe(flows = flow_missing, modeldata = md),
    "EpiSewer will impute missing flow data for 5 unsampled dates"
  )
  expect_equal(
    length(md_flow$flow), md_flow$T + md_flow$h
  )
  expect_true(
    all(!is.na(md_flow$flow))
  )
}
)

test_that("flows_observe() imputes flow values missing for future dates", {
  # take 5 unsampled dates
  flow_missing <- ww_data_SARS_CoV_2_Zurich$flows[date <= "2022-04-07",]
  expect_message(
    md_flow <- flows_observe(flows = flow_missing, modeldata = md),
    "EpiSewer will impute missing flow data for 7 forecasting dates"
  )
  expect_equal(
    length(md_flow$flow), md_flow$T + md_flow$h
  )
  expect_true(
    all(!is.na(md_flow$flow))
  )
}
)

test_that("flows_observe() imputes flow values missing for unsampled and future dates", {
  # take 5 unsampled dates
  missing_dates <- 1 + ww_data_SARS_CoV_2_Zurich$measurements$date[c(1, 5, 7, 8, 12)]
  flow_missing <- ww_data_SARS_CoV_2_Zurich$flows[!(date %in% missing_dates) & date <= "2022-04-07",]
  expect_message(
    md_flow <- flows_observe(flows = flow_missing, modeldata = md),
    "EpiSewer will impute missing flow data for 5 unsampled and 7 forecasting dates"
  )
  expect_equal(
    length(md_flow$flow), md_flow$T + md_flow$h
  )
  expect_true(
    all(!is.na(md_flow$flow))
  )
}
)


