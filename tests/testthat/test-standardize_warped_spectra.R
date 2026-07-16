test_that(
  "standardise_warped_spectra aligns orthophosphate peak",
  {


    spectra <- load_test_spectra_matrix(
      n = 5
    )


    ppm <- as.numeric(
      colnames(spectra)
    )


    metadata <- data.frame(

      Spectra.ID = rownames(spectra),

      Landuse = factor(
        c(
          "Forest",
          "Forest",
          "Pasture",
          "Pasture",
          "Forest"
        )
      )

    )


    result <- NMRDivR::standardise_warped_spectra(

      ppm = ppm,

      spectra_warped = spectra,

      reference = rownames(spectra)[1],

      metadata = metadata,

      treatment = "Landuse"

    )


    expect_s3_class(
      result,
      "spectra_standardised"
    )


    expect_named(
      result,
      c(
        "standardized.spectra",
        "summary.stats",
        "distance.matrices"
      )
    )


    expect_equal(

      dim(result$standardized.spectra),

      dim(spectra)

    )


    expect_true(
      all(
        is.finite(
          result$standardized.spectra
        )
      )
    )


    expect_true(
      "Warped" %in%
        result$summary.stats$Dataset
    )


    expect_true(
      "Standardised" %in%
        result$summary.stats$Dataset
    )


  })
