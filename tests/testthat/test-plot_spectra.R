test_that("plot_spectra generates plotly objects", {

  ppm <- seq(10, -10, length.out = 100)

  spectra <- matrix(
    abs(rnorm(300)),
    nrow = 3,
    ncol = 100
  )

  colnames(spectra) <- ppm

  rownames(spectra) <- c(
    "Sample1",
    "Sample2",
    "Sample3"
  )


  result <- plot_spectra(
    spectra.matrix.raw = spectra
  )


  # Returns named list
  expect_type(result, "list")

  expect_length(result, 3)


  # Check each element is a plotly htmlwidget
  expect_true(
    all(
      purrr::map_lgl(
        result,
        function(x) {
          inherits(x, "plotly")
        }
      )
    )
  )

})


test_that("plot_spectra handles raw warped and standardised spectra", {

  ppm <- seq(10, -10, length.out = 100)


  raw <- matrix(
    abs(rnorm(300)),
    nrow = 3,
    ncol = 100
  )

  warped <- raw + matrix(
    rnorm(300, sd = 0.01),
    nrow = 3
  )

  standardised <- warped / max(warped)


  colnames(raw) <- ppm
  colnames(warped) <- ppm
  colnames(standardised) <- ppm


  rownames(raw) <- c(
    "Sample1",
    "Sample2",
    "Sample3"
  )

  rownames(warped) <- rownames(raw)
  rownames(standardised) <- rownames(raw)


  plots <- plot_spectra(
    spectra.matrix.raw = raw,
    spectra.matrix.warped = warped,
    spectra.matrix.standardised = standardised
  )


  expect_type(
    plots,
    "list"
  )

  expect_length(
    plots,
    3
  )


  expect_true(
    all(
      purrr::map_lgl(
        plots,
        function(x) {
          inherits(x, "plotly")
        }
      )
    )
  )


})


test_that("plot_spectra rejects empty input", {


  expect_error(
    plot_spectra(),
    "Provide at least one spectra matrix"
  )


})


test_that("plot_spectra preserves sample names", {


  ppm <- seq(5, -5, length.out = 50)

  spectra <- matrix(
    runif(100),
    nrow = 2,
    ncol = 50
  )

  colnames(spectra) <- ppm

  rownames(spectra) <- c(
    "N10Y23",
    "N11Y23"
  )


  plots <- plot_spectra(
    spectra.matrix.raw = spectra
  )


  expect_true(
    all(
      c("N10Y23", "N11Y23") %in% names(plots)
    )
  )


})
