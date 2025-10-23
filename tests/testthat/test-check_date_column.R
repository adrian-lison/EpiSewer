library(data.table)

test_that("check_date_column() throws error when date column does not exist", {
  dt <- data.table(
    other_col = c(1, 2, 3),
    another_col = c("A", "B", "C")
  )

  expect_error(
    check_date_column(dt, "date"),
    "The date column does not exist in the data"
  )

  expect_error(
    check_date_column(dt, "date"),
    "Expected a column named 'date', but the data contains: other_col, another_col"
  )

  expect_error(
    check_date_column(dt, "my_date_column"),
    "Please ensure the column `my_date_column` is correctly specified"
  )
})

test_that("check_date_column() accepts valid Date objects", {
  dt <- data.table(
    date = as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")),
    dummy = c(1, 2, 3)
    )
  result <- check_date_column(dt, "date")

  expect_s3_class(result[["date"]], "Date")
  expect_equal(length(result[["date"]]), 3)
  expect_equal(result[["date"]], as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")))
})

test_that("check_date_column() converts POSIXct to Date with message", {
  dt <- data.table(date = as.POSIXct(c("2022-01-01 12:00:00", "2022-01-02 12:00:00", "2022-01-03 12:00:00"), tz = "UTC"))

  expect_message(
    result <- check_date_column(dt, "date"),
    "The date column `date` was of type POSIXct/POSIXlt and has been converted to Date"
  )
  expect_s3_class(result[["date"]], "Date")
  expect_equal(result[["date"]], as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")))
})

test_that("check_date_column() converts POSIXlt to Date with message", {
  # Create POSIXct first, then convert to POSIXlt to avoid data.table warning
  dt <- as.data.table(data.frame(date = as.POSIXlt(as.POSIXct(
    c("2022-01-01 12:00:00", "2022-01-02 12:00:00", "2022-01-03 12:00:00"),
    tz = "UTC"
    ))))

  expect_message(
    result <- check_date_column(dt, "date"),
    "The date column `date` was of type POSIXct/POSIXlt and has been converted to Date"
  )

  expect_s3_class(result[["date"]], "Date")
  expect_equal(result[["date"]], as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")))
})

test_that("check_date_column() converts character dates", {
  dt <- data.table(date = c("2022-01-01", "2022-01-02", "2022-01-03"))
  result <- check_date_column(dt, "date")
  expect_s3_class(result[["date"]], "Date")
  expect_equal(result[["date"]], as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")))
})

test_that("check_date_column() converts factor dates", {
  dt <- data.table(date = as.factor(c("2022-01-01", "2022-01-02", "2022-01-03")))
  result <- check_date_column(dt, "date")
  expect_s3_class(result[["date"]], "Date")
  expect_equal(result[["date"]], as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")))
})

test_that("check_date_column() throws error for unparseable character strings", {
  dt <- data.table(date = c("not a date", "also not a date", "2022-01-03"))

  expect_warning(expect_error(
    check_date_column(dt, "date"),
    "Please check your date data for invalid values"
  ))
})

test_that("check_date_column() throws error for numeric types", {
  dt <- data.table(date = c(1, 2, 3))

  expect_error(
    check_date_column(dt, "date"),
    paste(
      "The date column `date` must be of type Date, POSIXct, POSIXlt,",
      "or a character string"
    )
  )

  expect_error(
    check_date_column(dt, "date"),
    "Found type: numeric"
  )
})

test_that("check_date_column() throws error for logical types", {
  dt <- data.table(date = c(TRUE, FALSE, TRUE))

  expect_error(
    check_date_column(dt, "date"),
    paste(
      "The date column `date` must be of type Date, POSIXct, POSIXlt,",
      "or a character string"
    )
  )

  expect_error(
    check_date_column(dt, "date"),
    "Found type: logical"
  )
})

test_that("check_date_column() throws error when conversion results in NA values", {
  dt <- data.table(date = as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")))
  dt[2, date := NA]

  expect_error(
    check_date_column(dt, "date"),
    "After conversion, 1 date value\\(s\\) in column `date` are NA or invalid"
  )
})

test_that("check_date_column() throws error when multiple NA values after conversion", {
  dt <- data.table(date = as.Date(c("2022-01-01", "2022-01-02", "2022-01-03", "2022-01-04")))
  dt[c(2, 4), date := NA]

  expect_error(
    check_date_column(dt, "date"),
    "After conversion, 2 date value\\(s\\) in column `date` are NA or invalid"
  )
})

test_that("check_date_column() uses custom column name in error messages", {
  dt <- data.table(date = c(1, 2, 3))

  expect_error(
    check_date_column(dt, "custom_date_col"),
    "The date column `custom_date_col` must be of type Date"
  )
})

test_that("check_date_column() preserves other columns in data.table", {
  dt <- data.table(
    date = as.Date(c("2022-01-01", "2022-01-02", "2022-01-03")),
    value = c(10, 20, 30),
    name = c("A", "B", "C")
  )

  result <- check_date_column(dt, "date")

  expect_equal(ncol(result), 3)
  expect_equal(result$value, c(10, 20, 30))
  expect_equal(result$name, c("A", "B", "C"))
})

test_that("check_date_column() handles POSIXct with timezone information", {
  dt <- data.table(date = as.POSIXct(c("2022-01-01", "2022-01-02"), tz = "UTC"))

  expect_message(
    result <- check_date_column(dt, "date"),
    "POSIXct/POSIXlt"
  )

  expect_s3_class(result[["date"]], "Date")
})
