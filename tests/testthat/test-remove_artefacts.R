test_that("remove_artefacts performs baseline correction", {

  spectra.matrix <- load_test_spectra_matrix()


  result <- remove_artefacts(
    spectra.matrix.raw = spectra.matrix,
    BaselineCorrect = TRUE,
    ZeroFill = TRUE,
    nperm = 99
  )


  expect_type(
    result,
    "list"
  )


  expect_named(
    result,
    c(
      "SpectraCorrected",
      "BaselineCorrection",
      "Procrustes",
      "ProcrustesTest",
      "ProcrustesPlot"
    )
  )


  # dimensions preserved
  expect_equal(
    dim(result$SpectraCorrected),
    dim(spectra.matrix)
  )


  # negative values removed
  expect_true(
    all(result$SpectraCorrected >= 0)
  )


})


test_that("remove_artefacts can zero negative values without baseline correction", {


  spectra.matrix <- load_test_spectra_matrix()


  expect_warning(

    result <- remove_artefacts(
      spectra.matrix.raw = spectra.matrix,
      BaselineCorrect = FALSE,
      ZeroFill = TRUE,
      nperm = 99
    ),

    "BaselineCorrection is advised"

  )


  expect_true(
    all(result$SpectraCorrected >= 0)
  )


})


test_that("remove_artefacts generates valid Procrustes diagnostics", {


  spectra.matrix <- load_test_spectra_matrix()


  result <- remove_artefacts(
    spectra.matrix,
    nperm = 99
  )


  expect_s3_class(
    result$Procrustes,
    "procrustes"
  )


  expect_s3_class(
    result$ProcrustesTest,
    "protest"
  )


  expect_s3_class(
    result$ProcrustesPlot,
    "ggplot"
  )


  expect_true(
    is.numeric(result$ProcrustesTest$t0)
  )


})
