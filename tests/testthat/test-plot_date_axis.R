date_range <- function(from, days) {
  as.Date(from) + c(0, days)
}

test_that("date_breaks_adaptive() uses monthly breaks for short ranges", {
  breaks <- date_breaks_adaptive()(date_range("2022-03-15", 240))
  expect_equal(breaks, seq(as.Date("2022-04-01"), by = "month", length.out = 8))
})

test_that("date_breaks_adaptive() widens breaks for long ranges", {
  f <- date_breaks_adaptive(max_breaks = 8)
  for (days in c(400, 800, 1500, 4000)) {
    limits <- date_range("2022-03-15", days)
    breaks <- f(limits)
    expect_lte(length(breaks), 8)
    expect_gte(length(breaks), 2)
    expect_true(all(breaks >= limits[1] & breaks <= limits[2]))
    expect_true(all(format(breaks, "%d") == "01"))
  }
  # long ranges switch to yearly breaks on January 1st
  breaks <- f(date_range("2015-06-01", 4000))
  expect_true(all(format(breaks, "%m-%d") == "01-01"))
})

test_that("date_breaks_adaptive() handles very short and missing limits", {
  f <- date_breaks_adaptive()
  breaks <- f(date_range("2022-03-15", 30))
  expect_gte(length(breaks), 2)
  expect_null(f(as.Date(c(NA, "2022-03-15"))))
})

test_that("date_labels_adaptive() only shows the year when it changes", {
  breaks <- as.Date(c("2022-11-01", "2022-12-01", "2023-01-01", "2023-02-01"))
  month <- format(breaks, "%b") # locale-dependent month abbreviations
  expect_equal(
    date_labels_adaptive(breaks),
    c(
      paste0(month[1], "\n2022"), month[2],
      paste0(month[3], "\n2023"), month[4]
    )
  )
  expect_equal(
    date_labels_adaptive(as.Date(c("2020-01-01", "2022-01-01", NA))),
    c("2020", "2022", NA)
  )
})
