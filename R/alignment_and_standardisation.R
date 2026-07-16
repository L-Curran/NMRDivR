#' Perform spectrum alignment using PepsNMR warping
#'
#' Wraps \code{PepsNMR::Warping} with defaults tuned for environmental spectra.
#' Evaluates structural distortion using Bray-Curtis distances and Procrustes analysis.
#'
#' @param spectra_matrix_raw Numeric matrix of unaligned spectra (samples × ppm).
#' @param normalization.type Normalisation method passed to PepsNMR::Warping.
#' @param reference.choice Reference selection strategy for warping.
#' @param nperm Number of permutations for Procrustes test.
#' @param returnWarpFunc Logical; return warping function.
#' @param returnReference Logical; return reference spectrum.
#' @param ... Additional arguments passed to PepsNMR::Warping.
#'
#' @return A list with warped spectra, diagnostics, and Procrustes results.
#'
#' @export
warp_me <- function(
    spectra_matrix_raw,
    normalization.type = "mean",
    reference.choice = "after",
    nperm = 999,
    returnWarpFunc = TRUE,
    returnReference = TRUE,
    ...
) {

  args <- list(...)

  ppm_grid <- colnames(spectra_matrix_raw)

  warp_output <- do.call(
    PepsNMR::Warping,
    c(
      list(
        Spectrum_data = spectra_matrix_raw,
        normalization.type = normalization.type,
        reference.choice = reference.choice,
        verbose = TRUE
      ),
      args
    )
  )

  spectra_warped <- if (is.list(warp_output) && !is.null(warp_output$Spectrum_data)) {
    warp_output$Spectrum_data
  } else {
    warp_output
  }

  bc_raw <- vegan::vegdist(spectra_matrix_raw, method = "bray")
  bc_warped <- vegan::vegdist(spectra_warped, method = "bray")

  pcoa_raw <- stats::cmdscale(bc_raw, k = 2)
  pcoa_warped <- stats::cmdscale(bc_warped, k = 2)

  proc <- vegan::procrustes(pcoa_raw, pcoa_warped, symmetric = TRUE)
  proc_test <- vegan::protest(pcoa_raw, pcoa_warped, permutations = nperm)

  df_proc <- data.frame(
    Raw_X = proc$X[, 1],
    Raw_Y = proc$X[, 2],
    Warped_X = proc$Yrot[, 1],
    Warped_Y = proc$Yrot[, 2]
  )

  df_long <- tidyr::pivot_longer(
    df_proc,
    cols = everything(),
    names_to = c("Dataset", ".value"),
    names_pattern = "(Raw|Warped)_(X|Y)"
  )

  pal <- c("#0072B2", "#D55E00")

  p <- ggplot2::ggplot(df_long, ggplot2::aes(X, Y, color = Dataset)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_segment(
      data = df_proc,
      ggplot2::aes(x = Raw_X, y = Raw_Y, xend = Warped_X, yend = Warped_Y),
      inherit.aes = FALSE,
      colour = "grey60",
      linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "Procrustes Analysis: Raw vs Warped Spectra",
      subtitle = paste0(
        "R2 = ", round(proc_test$t0, 3),
        ", p = ", signif(proc_test$signif, 3)
      )
    )

  structure(
    list(
      SpectraRaw = spectra_matrix_raw,
      SpectraWarped = spectra_warped,
      ppm_grid = ppm_grid,
      Warping = warp_output,
      Procrustes = proc,
      ProcrustesTest = proc_test,
      ProcrustesPlot = p
    ),
    class = "spectra_warped"
  )
}


#' Standardise warped spectra using chemical shift calibration
#'
#' Aligns orthophosphate peak to a target ppm and evaluates structural change.
#'
#' @param ppm Numeric ppm vector.
#' @param spectra_warped Matrix of warped spectra.
#' @param reference Reference spectrum name.
#' @param ortho_peak_range Numeric length-2 ppm range.
#' @param target_ppm Target ppm position.
#' @param metadata Sample metadata (must include Spectra.ID).
#' @param treatment Grouping variable in metadata.
#' @param treatment_colours Optional colour vector.
#' @param nperm PERMANOVA permutations.
#'
#' @return List with corrected spectra, stats, models, and plots.
#'
#' @export
standardise_warped_spectra <- function(
    ppm,
    spectra_warped,
    reference,
    ortho_peak_range = c(4.5, 6.5),
    target_ppm = 6.0,
    metadata,
    treatment,
    treatment_colours = NULL,
    nperm = 999
) {

  metadata <- metadata[!is.na(metadata[[treatment]]), , drop = FALSE]
  metadata <- metadata[metadata$Spectra.ID %in% rownames(spectra_warped), , drop = FALSE]

  common <- intersect(rownames(spectra_warped), metadata$Spectra.ID)

  spectra_warped <- spectra_warped[common, , drop = FALSE]
  metadata <- metadata[match(common, metadata$Spectra.ID), , drop = FALSE]

  peak_region <- ppm >= ortho_peak_range[1] & ppm <= ortho_peak_range[2]

  corrected_spectra <- t(apply(spectra_warped, 1, function(spectrum) {
    peak_ppm <- ppm[peak_region][which.max(spectrum[peak_region])]
    shift <- peak_ppm - target_ppm
    stats::approx(ppm - shift, spectrum, xout = ppm, rule = 2)$y
  }))

  colnames(corrected_spectra) <- ppm
  rownames(corrected_spectra) <- rownames(spectra_warped)

  warped_dist <- vegan::vegdist(spectra_warped, "bray")
  std_dist <- vegan::vegdist(corrected_spectra, "bray")

  summary_stats <- data.frame(
    Dataset = c("Warped", "Standardised"),
    Mean = c(mean(warped_dist), mean(std_dist)),
    Median = c(stats::median(warped_dist), stats::median(std_dist)),
    Variance = c(stats::var(warped_dist), stats::var(std_dist))
  )

  structure(
    list(
      standardized.spectra = corrected_spectra,
      summary.stats = summary_stats,
      distance.matrices = list(warped = warped_dist, standardized = std_dist)
    ),
    class = "spectra_standardised"
  )
}

#' Plot raw, warped, and standardised spectra
#'
#' Generates interactive Plotly spectra comparisons.
#'
#' @param spectra.matrix.raw Raw spectra matrix.
#' @param spectra.matrix.warped Warped spectra matrix.
#' @param spectra.matrix.standardised Standardised spectra matrix.
#' @param warping.reference Optional reference spectrum.
#' @param sample.names Optional subset of samples.
#'
#' @return Named list of Plotly objects.
#'
#' @export
plot_spectra <- function(
    spectra.matrix.raw = NULL,
    spectra.matrix.warped = NULL,
    spectra.matrix.standardised = NULL,
    warping.reference = NULL,
    sample.names = NULL
) {

  if (is.null(spectra.matrix.raw) &&
      is.null(spectra.matrix.warped) &&
      is.null(spectra.matrix.standardised)) {
    stop("Provide at least one spectra matrix.")
  }

  build_list <- function(mat, name) {
    if (is.null(mat)) return(NULL)
    lapply(seq_len(nrow(mat)), function(i) {
      list(
        name = rownames(mat)[i],
        ppm = as.numeric(colnames(mat)),
        intensity = as.numeric(mat[i, ]),
        type = name
      )
    })
  }

  all_spectra <- c(
    build_list(spectra.matrix.raw, "Raw"),
    build_list(spectra.matrix.warped, "Warped"),
    build_list(spectra.matrix.standardised, "standardised")
  )

  sample_names <- unique(vapply(all_spectra, `[[`, "", "name"))

  plots <- list()

  for (s in sample_names) {

    sample_data <- Filter(function(x) x$name == s, all_spectra)

    p <- plotly::plot_ly()

    for (sp in sample_data) {
      p <- plotly::add_trace(
        p,
        x = sp$ppm,
        y = sp$intensity,
        type = "scatter",
        mode = "lines",
        name = sp$type
      )
    }

    p <- plotly::layout(
      p,
      title = paste0("Spectra: ", s),
      xaxis = list(title = "PPM", autorange = "reversed"),
      yaxis = list(title = "Intensity")
    )

    plots[[s]] <- p
  }

  plots
}
