test_that("warp_me returns valid warped spectra object", {


  spectra <- load_test_spectra_matrix(
    n = 5
  )


  result <- NMRDivR::warp_me(
    spectra,
    nperm = 99
  )


  expect_s3_class(
    result,
    "spectra_warped"
  )


  expect_true(
    is.list(result)
  )


  expect_named(
    result,
    c(
      "SpectraRaw",
      "SpectraWarped",
      "ppm_grid",
      "Warping",
      "Procrustes",
      "ProcrustesTest",
      "ProcrustesPlot"
    )
  )


  expect_equal(
    dim(result$SpectraRaw),
    dim(result$SpectraWarped)
  )


  expect_equal(
    colnames(result$SpectraRaw),
    colnames(result$SpectraWarped)
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

})

test_that("warp_me rejects invalid spectra input", {


  expect_error(
    NMRDivR::warp_me(
      matrix(
        letters[1:10],
        nrow=2
      )
    )
  )


})
