test_that("flatten_one_level removes one level and keeps inner names", {
  x <- list(a = list(x = 1, y = "b"), b = 3, c = list(z = list(1, 2)))
  expect_identical(
    flatten_one_level(x),
    list(x = 1, y = "b", b = 3, z = list(1, 2))
  )
  expect_identical(flatten_one_level(list()), list())
})

test_that("flatten_one_level keeps data.frames intact", {
  dt <- data.table::data.table(a = 1:2)
  out <- flatten_one_level(list(m = list(d = dt), e = dt))
  expect_s3_class(out$d, "data.table")
  expect_s3_class(out$e, "data.table")
})

test_that("assert_valid_summary rejects failed or missing summaries", {
  expect_error(assert_valid_summary(list()), "no valid summary")
  failed <- list(summary = structure("boom", class = "try-error"))
  expect_error(assert_valid_summary(failed), "boom")
  expect_true(assert_valid_summary(list(summary = list(R = 1))))
  expect_error(combine_summaries(list(failed), "R"), "no valid summary")
})
