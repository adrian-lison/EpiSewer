# One-off script to record modeldata characterization fixtures.
#
# Run manually from the package root on a KNOWN-GOOD commit:
#   Rscript data-raw/record_modeldata_fixtures.R
#
# For each configuration in the battery (tests/testthat/helper-modeldata_configs.R),
# this builds the EpiSewer job without fitting and stores the stan data, inits
# and the externally-consumed metainfo fields. The comparison tests in
# test-modeldata_characterization.R rebuild the jobs and check equality.
#
# The recorded fixtures live in tests/testthat/fixtures/modeldata/, which is
# gitignored: they are a local migration gate, and the comparison tests skip
# when fixtures are absent (e.g. on CI).

suppressMessages(devtools::load_all(quiet = TRUE))
source("tests/testthat/helper-modeldata_configs.R")

out_dir <- file.path("tests", "testthat", "fixtures", "modeldata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

battery <- modeldata_config_battery()

for (name in names(battery)) {
  cat("Recording", name, "... ")
  fixture <- suppressWarnings(suppressMessages(
    modeldata_fixture_job(battery[[name]]())
  ))
  saveRDS(
    fixture,
    file.path(out_dir, paste0(name, ".rds")),
    version = 2, compress = "xz"
  )
  cat("ok (", length(fixture$data), "data vars )\n")
}

cat("Done. Recorded", length(battery), "fixtures in", out_dir, "\n")
