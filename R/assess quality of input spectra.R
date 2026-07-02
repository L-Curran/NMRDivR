#' Perform baseline correction and negative value zeroing
#'
#' This function applies optional baseline correction and negative value
#' removal to a spectral matrix. It optionally uses PepsNMR preprocessing
#' methods and evaluates structural preservation using Bray-Curtis distances
#' and Procrustes analysis.
#'
#' @param spectra.matrix.raw Numeric matrix of raw spectra
#' (samples × ppm).
#'
#' @param BaselineCorrect Logical; whether to apply baseline correction
#' using `PepsNMR::BaselineCorrection`.
#'
#' @param ZeroFill Logical; whether to replace negative values with zero
#' using `PepsNMR::NegativeValuesZeroing`.
#'
#' @param nperm Number of permutations for Procrustes significance test.
#'
#' @param ... Additional arguments passed to PepsNMR functions.
#'
#' @return A list containing:
#' \item{SpectraCorrected}{Processed spectral matrix.}
#' \item{BaselineCorrection}{Baseline correction output (if applied).}
#' \item{Procrustes}{Procrustes analysis object.}
#' \item{ProcrustesTest}{Permutation test results.}
#' \item{ProcrustesPlot}{ggplot diagnostic plot.}
#'
#' @export
remove_artefacts <- function(spectra.matrix.raw,
                             BaselineCorrect = TRUE,
                             ZeroFill = TRUE,
                             nperm = 999,
                             ...) {

  args <- list(...)
  spectra_corrected <- spectra.matrix.raw
  baseline_output <- NULL

  #=========================================================
  # BASELINE CORRECTION
  #=========================================================
  if (BaselineCorrect) {

    baseline_output <- do.call(
      PepsNMR::BaselineCorrection,
      c(
        list(
          Spectrum_data = spectra_corrected,
          ptw.bc = TRUE,
          lambda.bc = 1e7,
          ppm.bc = TRUE,
          returnBaseline = TRUE,
          verbose = TRUE
        ),
        args
      )
    )

    spectra_corrected <- baseline_output$Spectrum_data
  }

  #=========================================================
  # ZEROING NEGATIVE VALUES
  #=========================================================
  if (ZeroFill) {

    if (!BaselineCorrect) {
      warning("BaselineCorrection is advised before ZeroFilling")
    }

    spectra_corrected <- do.call(
      PepsNMR::NegativeValuesZeroing,
      c(
        list(Spectrum_data = spectra_corrected, verbose = TRUE),
        args
      )
    )
  }

  #=========================================================
  # DISTANCES (BRAY-CURTIS)
  #=========================================================
  bc_raw <- vegan::vegdist(spectra.matrix.raw, method = "bray")
  bc_corr <- vegan::vegdist(spectra_corrected, method = "bray")

  #=========================================================
  # PCoA (cmdscale)
  #=========================================================
  pcoa_raw <- stats::cmdscale(bc_raw, k = 2)
  pcoa_corr <- stats::cmdscale(bc_corr, k = 2)

  #=========================================================
  # PROCRUSTES ANALYSIS
  #=========================================================
  proc <- vegan::procrustes(pcoa_raw, pcoa_corr, symmetric = TRUE)
  proc_test <- vegan::protest(pcoa_raw, pcoa_corr, permutations = nperm)

  #=========================================================
  # DATA FRAME (more efficient construction)
  #=========================================================
  df_proc <- data.frame(
    Raw_X = proc$X[, 1],
    Raw_Y = proc$X[, 2],
    Corrected_X = proc$Yrot[, 1],
    Corrected_Y = proc$Yrot[, 2],
    Spectrum = rownames(proc$X)
  )

  df_long <- tidyr::pivot_longer(
    df_proc,
    cols = c(Raw_X, Raw_Y, Corrected_X, Corrected_Y),
    names_to = c("Dataset", ".value"),
    names_pattern = "(Raw|Corrected)_(X|Y)"
  )

  #=========================================================
  # PLOT (fully namespaced + simplified)
  #=========================================================
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = X, y = Y, color = Dataset)) +

    ggplot2::stat_ellipse(
      ggplot2::aes(group = Dataset, fill = Dataset),
      geom = "polygon",
      alpha = 0.05
    ) +

    ggplot2::geom_point(size = 3) +

    ggplot2::geom_segment(
      data = df_proc,
      ggplot2::aes(x = Raw_X, y = Raw_Y,
                   xend = Corrected_X, yend = Corrected_Y),
      inherit.aes = FALSE,
      color = "grey50",
      linetype = "dashed"
    ) +

    ggplot2::scale_color_manual(values = c(Raw = "#0072B2", Corrected = "#D55E00")) +
    ggplot2::scale_fill_manual(values = c(Raw = "#0072B2", Corrected = "#D55E00")) +

    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      axis.title = ggplot2::element_blank()
    ) +

    ggplot2::labs(
      title = "Procrustes Analysis: Raw vs Corrected Spectra",
      subtitle = paste0(
        "R² = ", round(proc_test$t0, 3),
        ", p = ", signif(proc_test$signif, 3)
      )
    )

  #=========================================================
  # RETURN
  #=========================================================
  list(
    SpectraCorrected = spectra_corrected,
    BaselineCorrection = baseline_output,
    Procrustes = proc,
    ProcrustesTest = proc_test,
    ProcrustesPlot = p
  )
}

#' Assess quality of manually identified peaks
#'
#' Evaluates manually annotated peaks for coverage and positional variance.
#' Filters unstable or poorly represented peaks and visualises peak
#' shift distributions.
#'
#' @param metadata Data frame containing peak annotations.
#'
#' @param sample.id Column name for sample identifiers.
#'
#' @param peak.id Column name for peak identifiers.
#'
#' @param shift Column name for chemical shift (ppm).
#'
#' @param coverage Minimum proportion of samples a peak must appear in.
#'
#' @param sd Multiplier for filtering peaks based on shift variance.
#'
#' @return A list containing filtered metadata, statistics, and plots.
#'
#' @export
assess_manual_peak.ids <- function(metadata,
                                   sample.id = "Spectra.ID",
                                   peak.id = "Peak.ID",
                                   shift = "Peak.Shift",
                                   coverage = 0.8,
                                   sd = 2) {

  metadata <- as.data.frame(metadata)

  total_samples <- length(unique(metadata[[sample.id]]))
  all_peaks <- unique(metadata[[peak.id]])

  # -----------------------------
  # Peak statistics (fast dplyr)
  # -----------------------------
  peak_stats <- dplyr::group_by(metadata, .data[[peak.id]]) |>
    dplyr::summarise(
      mean_shift = mean(.data[[shift]], na.rm = TRUE),
      var_shift  = stats::var(.data[[shift]], na.rm = TRUE),
      sd_shift   = stats::sd(.data[[shift]], na.rm = TRUE),
      n_samples  = dplyr::n_distinct(.data[[sample.id]]),
      coverage   = n_samples / total_samples,
      .groups = "drop"
    )

  sd_threshold <- sd * stats::median(peak_stats$sd_shift, na.rm = TRUE)

  # -----------------------------
  # Filtering
  # -----------------------------
  valid_peaks <- peak_stats[[peak.id]][
    peak_stats$coverage >= coverage &
      peak_stats$sd_shift <= sd_threshold
  ]

  filtered_metadata <- metadata[metadata[[peak.id]] %in% valid_peaks, , drop = FALSE]
  filtered_stats <- peak_stats[peak_stats[[peak.id]] %in% valid_peaks, , drop = FALSE]

  filtered_stats <- filtered_stats[order(filtered_stats$sd_shift), , drop = FALSE]

  filtered_metadata$order_peak <- factor(
    filtered_metadata[[peak.id]],
    levels = filtered_stats[[peak.id]]
  )

  # -----------------------------
  # Plot (fixed geom_vline usage)
  # -----------------------------
  plot <- ggplot2::ggplot(filtered_metadata,
                          ggplot2::aes(x = .data[[shift]])) +

    ggplot2::geom_density(fill = "steelblue", alpha = 0.4) +

    ggplot2::facet_wrap(~ order_peak, scales = "free") +

    ggplot2::geom_vline(
      data = filtered_stats,
      ggplot2::aes(xintercept = mean_shift - 2 * sd_shift),
      linetype = "dashed",
      colour = "red"
    ) +

    ggplot2::geom_vline(
      data = filtered_stats,
      ggplot2::aes(xintercept = mean_shift + 2 * sd_shift),
      linetype = "dashed",
      colour = "red"
    ) +

    ggplot2::geom_vline(
      data = filtered_stats,
      ggplot2::aes(xintercept = mean_shift),
      colour = "black"
    ) +

    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Peak Shift distributions (SD-normalised filtering)",
      x = "Peak Shift",
      y = "Density"
    )

  print(plot)

  cat(
    "PEAK FILTER SUMMARY\n",
    "-------------------\n",
    "Original Peak IDs:", length(all_peaks), "\n",
    "Filtered Peak IDs:", length(valid_peaks), "\n",
    "SD threshold:", round(sd_threshold, 4), "\n"
  )

  list(
    input_metadata = metadata,
    peak_statistics = peak_stats,
    filtered_metadata = filtered_metadata,
    filtered_statistics = filtered_stats,
    sd_threshold = sd_threshold,
    plot = plot
  )
}

#' Inspect raw and corrected spectra
#'
#' Generates an interactive Plotly overlay of raw and corrected spectra
#' to visually assess baseline artefacts and correction performance.
#'
#' @param spectra.matrix.raw Raw spectral matrix (samples × ppm).
#'
#' @param spectra.matrix.corrected Corrected spectral matrix.
#'
#' @return A Plotly object displaying spectral overlays.
#'
#' @export
inspect_raw_spectra <- function(spectra.matrix.raw = NULL,
                                spectra.matrix.corrected = NULL) {

  if (is.null(spectra.matrix.raw) &&
      is.null(spectra.matrix.corrected)) {
    stop("At least one matrix must be supplied.")
  }

  get_palette <- function(n) {
    cb <- c("#0072B2", "#D55E00", "#009E73",
            "#CC79A7", "#56B4E9", "#E69F00",
            "#F0E442")

    if (n <= length(cb)) cb[seq_len(n)] else viridisLite::viridis(n)
  }

  sample_names <- unique(c(
    rownames(spectra.matrix.raw),
    rownames(spectra.matrix.corrected)
  ))

  sample_cols <- setNames(get_palette(length(sample_names)), sample_names)

  p <- plotly::plot_ly()

  # ---------------- raw ----------------
  if (!is.null(spectra.matrix.raw)) {
    ppm <- as.numeric(colnames(spectra.matrix.raw))

    for (i in seq_len(nrow(spectra.matrix.raw))) {
      sample <- rownames(spectra.matrix.raw)[i]

      p <- plotly::add_trace(
        p,
        x = ppm,
        y = spectra.matrix.raw[i, ],
        type = "scattergl",
        mode = "lines",
        name = sample,
        legendgroup = sample,
        showlegend = TRUE,
        line = list(color = sample_cols[sample], dash = "dot", width = 1)
      )
    }
  }

  # ---------------- corrected ----------------
  if (!is.null(spectra.matrix.corrected)) {
    ppm <- as.numeric(colnames(spectra.matrix.corrected))

    for (i in seq_len(nrow(spectra.matrix.corrected))) {
      sample <- rownames(spectra.matrix.corrected)[i]

      p <- plotly::add_trace(
        p,
        x = ppm,
        y = spectra.matrix.corrected[i, ],
        type = "scattergl",
        mode = "lines",
        name = sample,
        legendgroup = sample,
        showlegend = FALSE,
        line = list(color = sample_cols[sample], dash = "solid", width = 2)
      )
    }
  }

  plotly::layout(
    p,
    title = "Raw (dashed) vs Corrected (solid) Spectra",
    xaxis = list(title = "ppm", autorange = "reversed"),
    yaxis = list(title = "Intensity"),
    legend = list(orientation = "v"),
    hovermode = "closest"
  )
}
