load(
  testthat::test_path(
    "test-data",
    "spectra_matrix_standardized.RData"
  )
)

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

ppm <- as.numeric(colnames(spectra.matrix.stand))

sample.names <- rownames(spectra.matrix.stand)

peaks <- detect_peaks(
  spectra.matrix = spectra.matrix.stand,
  ppm = ppm,
  sample.names = sample.names
)


############################################################
# detect_peaks()
############################################################


test_that(
  "detect_peaks validates input types",
  {

    ppm <- as.numeric(colnames(spectra.matrix.stand))
    sample.names <- rownames(spectra.matrix.stand)


    # spectra.matrix must be a matrix
    expect_error(
      detect_peaks(
        spectra.matrix = as.data.frame(spectra.matrix.stand),
        ppm = ppm,
        sample.names = sample.names
      ),
      "must be a matrix"
    )


    # ppm must be numeric
    expect_error(
      detect_peaks(
        spectra.matrix = spectra.matrix.stand,
        ppm = as.character(ppm),
        sample.names = sample.names
      ),
      "must be numeric"
    )


    # ppm length must match number of columns
    expect_error(
      detect_peaks(
        spectra.matrix = spectra.matrix.stand,
        ppm = ppm[-1],
        sample.names = sample.names
      ),
      "must equal"
    )


    # sample.names length must match rows
    expect_error(
      detect_peaks(
        spectra.matrix = spectra.matrix.stand,
        ppm = ppm,
        sample.names = sample.names[-1]
      ),
      "must equal"
    )

  }
)


test_that(
  "detect_peaks returns expected peak metadata",
  {

    ppm <- as.numeric(colnames(spectra.matrix.stand))

    sample.names <- rownames(spectra.matrix.stand)


    peaks <- detect_peaks(
      spectra.matrix = spectra.matrix.stand,
      ppm = ppm,
      sample.names = sample.names
    )


    expect_s3_class(
      peaks,
      "data.frame"
    )


    expect_true(
      nrow(peaks) > 0
    )


    expect_true(
      all(
        c(
          "Peak.ID",
          "Peak.Shift",
          "Peak.AUC",
          "Peak.SNR",
          "Peak.Scale",
          "Spectra.ID",
          "Peak.Width",
          "Peak.Shift.SD"
        ) %in% colnames(peaks)
      )
    )


    expect_true(
      all(peaks$Peak.Width > 0)
    )


    expect_true(
      all(peaks$Peak.Shift.SD > 0)
    )


    expect_true(
      all(peaks$Spectra.ID %in% sample.names)
    )

  }
)



############################################################
# peak_assessment()
############################################################



test_that(
  "peak_assessment validates required columns",
  {


    bad_df <- data.frame(
      Peak.SNR = 1:5
    )


    expect_error(
      peak_assessment(bad_df),
      "must contain columns"
    )


  })




test_that(
  "peak_assessment filters low quality peaks",
  {


    result <- peak_assessment(
      peaks,
      min_signal = 10
    )


    expect_true(
      is.list(result)
    )


    expect_true(
      all(
        result$filtered.peaks$Peak.SNR >= 10
      )
    )


    expect_true(
      nrow(result$filtered.peaks) <
        nrow(peaks)
    )


  })




test_that(
  "peak_assessment creates diagnostic plots",
  {


    result <- peak_assessment(
      peaks
    )


    expect_true(
      inherits(
        result$SNR.density.plot,
        "ggplot"
      )
    )


    expect_true(
      is.data.frame(
        result$peak.outliers
      )
    )


  })




test_that(
  "peak_assessment joins metadata correctly",
  {



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

    result <- peak_assessment(
      peaks,
      spectra_metadata = metadata,
      treatment = "Landuse"
    )



    expect_true(
      "Landuse" %in%
        colnames(result$filtered.peaks)
    )


    expect_true(
      all(
        !is.na(result$filtered.peaks$Landuse)
      )
    )


  })





############################################################
# plot_detected_peaks()
############################################################



test_that(
  "plot_detected_peaks validates dataframe",
  {


    bad_df <- data.frame(
      Sample="A"
    )


    expect_error(
      plot_detected_peaks(bad_df),
      "must contain"
    )


  })




test_that(
  "plot_detected_peaks generates plotly object",
  {


    peaks <- detect_peaks(
      spectra.matrix = spectra.matrix.stand,
      ppm = ppm,
      sample.names = sample.names
    )


    result <- plot_detected_peaks(
      peaks
    )



    expect_true(
      inherits(
        result,
        "plotly"
      )
    )


  })





test_that(
  "plot_detected_peaks handles spectra overlay",
  {

    p <- plot_detected_peaks(
      peaks,
      spectra_matrix = spectra.matrix.stand
    )

    expect_s3_class(
      p,
      "plotly"
    )

  }
)



test_that(
  "plot_detected_peaks handles metadata colouring",
  {


    peaks <- plot_detected_peaks(
      peaks,
      spectra_matrix = spectra.matrix.stand
    )



    rownames(metadata) <- metadata$Spectra.ID



    result <- plot_detected_peaks(
      peaks,
      metadata = metadata,
      treatment = "Landuse"
    )



    expect_true(
      inherits(
        result,
        "plotly"
      )
    )


  })
