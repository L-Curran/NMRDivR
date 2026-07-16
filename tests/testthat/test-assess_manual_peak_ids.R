test_that("manual peak metadata contains required columns", {


  metadata <- read.csv(
    testthat::test_path(
      "test-data",
      "metadata.csv"
    )
  )


  expect_true(
    all(
      c(
        "Spectra.ID",
        "Peak.ID",
        "Peak.Shift"
      ) %in% names(metadata)
    )
  )


})


test_that("assess_manual_peak.ids filters unstable peaks", {


  metadata <- read.csv(
    testthat::test_path(
      "test-data",
      "metadata.csv"
    )
  )


  result <- assess_manual_peak.ids(
    metadata,
    coverage = 0.5
  )


  expect_type(
    result,
    "list"
  )


  expect_named(
    result,
    c(
      "input_metadata",
      "peak_statistics",
      "filtered_metadata",
      "filtered_statistics",
      "sd_threshold",
      "plot"
    )
  )


  expect_true(
    nrow(result$filtered_metadata) > 0
  )


  expect_true(
    nrow(result$filtered_statistics) > 0
  )


  expect_true(
    is.numeric(result$sd_threshold)
  )


})


test_that("filtered peaks meet coverage threshold", {


  metadata <- read.csv(
    testthat::test_path(
      "test-data",
      "metadata.csv"
    )
  )


  coverage_threshold <- 0.75


  result <- assess_manual_peak.ids(
    metadata,
    coverage = coverage_threshold
  )


  expect_true(
    all(
      result$filtered_statistics$coverage >= coverage_threshold
    )
  )


})
