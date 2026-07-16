test_that("sample_sanity keeps matching samples", {

  spectra <- matrix(
    runif(30),
    nrow = 3
  )

  rownames(spectra) <- c(
    "A",
    "B",
    "C"
  )


  metadata <- data.frame(
    Spectra.ID = c(
      "A",
      "B",
      "C"
    ),
    Treatment = c(
      "x",
      "y",
      "z"
    )
  )


  result <- sample_sanity(
    spectra,
    metadata
  )


  expect_equal(
    dim(result$spectra.matrix),
    dim(spectra)
  )


  expect_equal(
    nrow(result$spectra.meta),
    3
  )

})


test_that("sample_sanity removes spectra missing metadata", {

  spectra <- matrix(
    1:12,
    nrow = 3
  )

  rownames(spectra) <- c(
    "A",
    "B",
    "C"
  )


  metadata <- data.frame(
    Spectra.ID = c(
      "A",
      "B"
    )
  )


  result <- suppressWarnings(
    sample_sanity(
      spectra,
      metadata
    )
  )


  expect_equal(
    rownames(result$spectra.matrix),
    c(
      "A",
      "B"
    )
  )


  expect_equal(
    result$spectra.meta$Spectra.ID,
    c(
      "A",
      "B"
    )
  )

})

test_that("sample_sanity removes metadata without spectra", {

  spectra <- matrix(
    1:8,
    nrow = 2
  )

  rownames(spectra) <- c(
    "A",
    "B"
  )


  metadata <- data.frame(
    Spectra.ID = c(
      "A",
      "B",
      "C"
    )
  )


  result <- suppressWarnings(
    sample_sanity(
      spectra,
      metadata
    )
  )


  expect_equal(
    nrow(result$spectra.meta),
    2
  )


  expect_false(
    "C" %in% result$spectra.meta$Spectra.ID
  )

})


test_that("sample_sanity returns expected structure", {

  spectra <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(
      c("A","B"),
      NULL
    )
  )


  metadata <- data.frame(
    Spectra.ID = c(
      "A",
      "B"
    )
  )


  result <- sample_sanity(
    spectra,
    metadata
  )


  expect_named(
    result,
    c(
      "spectra.matrix",
      "spectra.meta"
    )
  )


  expect_true(
    is.matrix(result$spectra.matrix)
  )


  expect_true(
    is.data.frame(result$spectra.meta)
  )

})


test_that("test metadata file exists and loads", {

  metadata_file <- testthat::test_path(
    "test-data",
    "metadata.csv"
  )

  expect_true(
    file.exists(metadata_file)
  )


  metadata <- read.csv(
    metadata_file
  )


  expect_true(
    nrow(metadata) > 0
  )


  expect_true(
    "Spectra.ID" %in% names(metadata)
  )

})
