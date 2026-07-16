test_that("load_NMRData1D loads Bruker spectra correctly", {

  test_dir <- testthat::test_path("test-data")

  spectra_folder <- file.path(
    test_dir,
    "spectra"
  )

  samples <- list.dirs(
    spectra_folder,
    recursive = FALSE,
    full.names = FALSE
  )

  expect_true(length(samples) > 0)

  result <- load_NMRData1D(
    spectra.folder = spectra_folder,
    samples = samples,
    format = "bruker"
  )

  # Object structure
  expect_type(result, "list")

  # Named according to samples
  expect_setequal(
    names(result),
    samples
  )

  # No missing spectra
  expect_false(
    any(vapply(result, is.null, logical(1)))
  )

  # Each object should be an NMRData1D object
  expect_true(
    all(
      vapply(
        result,
        function(x) inherits(x, "NMRData1D"),
        logical(1)
      )
    )
  )

})


test_that("load_NMRData1D handles missing spectra", {

  test_dir <- testthat::test_path("test-data")

  spectra_folder <- file.path(
    test_dir,
    "spectra"
  )


  real_sample <- basename(
    list.dirs(
      spectra_folder,
      recursive = FALSE,
      full.names = TRUE
    )[1]
  )


  result <- suppressWarnings(
    load_NMRData1D(
      spectra.folder = spectra_folder,
      samples = c(
        real_sample,
        "missing_sample"
      ),
      format = "bruker"
    )
  )


  expect_true(
    "missing_sample" %in% names(result)
  )


  expect_null(
    result[["missing_sample"]]
  )


  expect_s4_class(
    result[[real_sample]],
    "NMRData1D"
  )

})
test_that("load_NMRData1D rejects invalid formats", {

  expect_error(

    load_NMRData1D(
      spectra.folder = "fake",
      samples = "sample",
      format = "invalid"
    ),

    regexp = "should be one of"
  )

})


test_that("JCAMP loading requires processing parameters", {

  test_dir <- tempfile()

  dir.create(
    file.path(
      test_dir,
      "sample",
      "2"
    ),
    recursive = TRUE
  )


  expect_error(

    load_NMRData1D(
      spectra.folder = test_dir,
      samples = "sample",
      format = "jcamp-dx"
    ),

    regexp = "JCAMP-DX requires"

  )

})
