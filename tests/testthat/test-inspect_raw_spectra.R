test_that("inspect_raw_spectra creates plotly object", {


  spectra.matrix <- load_test_spectra_matrix()


  result <- inspect_raw_spectra(
    spectra.matrix.raw = spectra.matrix
  )


  expect_s3_class(
    result,
    "plotly"
  )


})


test_that("inspect_raw_spectra accepts raw and corrected spectra", {


  spectra.matrix <- load_test_spectra_matrix()


  corrected <- spectra.matrix

  corrected[corrected < 0] <- 0


  result <- inspect_raw_spectra(
    spectra.matrix.raw = spectra.matrix,
    spectra.matrix.corrected = corrected
  )


  expect_s3_class(
    result,
    "plotly"
  )


})

test_that("inspect_raw_spectra requires input", {


  expect_error(

    inspect_raw_spectra(),

    "At least one matrix must be supplied"

  )


})
