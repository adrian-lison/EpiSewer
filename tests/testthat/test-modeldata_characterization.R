# Characterization tests for modeldata construction.
#
# Each configuration in the battery (helper-modeldata_configs.R) is rebuilt
# and compared against fixtures recorded on a known-good commit via
# fixtures/record_modeldata_fixtures.R. Comparison is on name-sorted content:
# element ORDER of job$data/init is allowed to change, values are not.

battery <- modeldata_config_battery()

for (config_name in names(battery)) {
  test_that(paste0("modeldata unchanged for config '", config_name, "'"), {
    fixture_file <- modeldata_fixture_path(config_name)
    skip_if_not(file.exists(fixture_file), "fixture not recorded")
    fixture <- readRDS(fixture_file)

    rebuilt <- suppressWarnings(suppressMessages(
      modeldata_fixture_job(battery[[config_name]]())
    ))

    expect_equal(names(rebuilt$data), names(fixture$data))
    expect_equal(rebuilt$data, fixture$data)
    expect_equal(names(rebuilt$init), names(fixture$init))
    expect_equal(rebuilt$init, fixture$init)
    expect_equal(rebuilt$metainfo, fixture$metainfo)
  })
}

# Pin the user-facing message texts of the message-heavy configurations.
# Message ORDER may change across refactors (single-pass vs multi-pass
# resolution); if only the order changed, updating the snapshot is expected.

test_that("messages for config 'calibrate_cases_direct' are stable", {
  expect_snapshot({
    invisible(modeldata_fixture_job(battery$calibrate_cases_direct()))
  })
})

test_that("messages for config 'default_sars' are stable", {
  expect_snapshot({
    invisible(modeldata_fixture_job(battery$default_sars()))
  })
})

test_that("messages for config 'partitions_likelihood' are stable", {
  expect_snapshot({
    invisible(modeldata_fixture_job(battery$partitions_likelihood()))
  })
})
