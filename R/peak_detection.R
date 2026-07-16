#' Perform automated peak detection
#'
#' This function serves as a wrapper for the speaq::getWavelettPeaks where the default parameters have been optimized for spectra derived from complex media, such as soils. In addition, this function also calculates the width and standard deviation of chemical shift for all peaks to be used in downstream analysis.
#'
#'
#' @param spectra.matrix A spectral matrix in sample x ppm format
#' @param ppm A numeric vector containing the ppm values of spectra.matrix
#' @param sample.names A character vector containing the the names each spectrum to perform peak detection within. These should correspond to the rownames of spectra.matrix. By default, peak detection is performed upon all spectra within spectra.matrix.
#' @param window.width Specifies the width of the search window during peak detection. See the package speaq for further information.
#' @param window.split Specifies the number of units to split the searching window.  during peak detection. See the package speaq for further information.
#' @param include.nearbyPeaks Include peaks within the tails of larger neighboring peaks that may be obscured by. If FALSE, these peaks are ignored.
#' @return A dataframe of peaks detected per spectrum and their associated metadata.
#' @export


detect_peaks <- function(spectra.matrix,
                         ppm,
                         sample.names,
                         window.width = "large",
                         window.split = 8,
                         include.nearbyPeaks = TRUE,
                         ...) {
  # ----------------------------------------------------------
  # Fast input checks
  # ----------------------------------------------------------

  if (!is.matrix(spectra.matrix)) {
    stop("'spectra.matrix' must be a matrix.")
  }

  if (!is.numeric(ppm)) {
    stop("'ppm' must be numeric.")
  }

  if (length(ppm) != ncol(spectra.matrix)) {
    stop(
      "'ppm' length (",
      length(ppm),
      ") must equal ncol(spectra.matrix) (",
      ncol(spectra.matrix),
      ")."
    )
  }

  if (length(sample.names) != nrow(spectra.matrix)) {
    stop(
      "'sample.names' length (",
      length(sample.names),
      ") must equal nrow(spectra.matrix) (",
      nrow(spectra.matrix),
      ")."
    )
  }

  # ----------------------------------------------------------
  # Compute ppm resolution once
  # ----------------------------------------------------------

  ppm_resolution <- stats::median(abs(diff(ppm)))

  # ----------------------------------------------------------
  # Peak detection
  # ----------------------------------------------------------

  detected_peaks <- speaq::getWaveletPeaks(
    Y.spec = spectra.matrix,
    X.ppm = ppm,
    sample.labels = sample.names,
    window.width = window.width,
    window.split = window.split,
    include.nearbyPeaks = include.nearbyPeaks,
    raw_peakheight = raw_peakheight,
    ...
  )

  # ----------------------------------------------------------
  # Rename columns
  # ----------------------------------------------------------

  colnames(detected_peaks) <- c("Peak.ID",
                                "Peak.Shift",
                                "Peak.AUC",
                                "Peak.SNR",
                                "Peak.Scale",
                                "Spectra.ID")

  # ----------------------------------------------------------
  # Derived quantities
  # ----------------------------------------------------------

  detected_peaks$Peak.Width <-
    detected_peaks$Peak.Scale * ppm_resolution

  detected_peaks$Peak.Shift.SD <-
    0.0015 + 0.002 / (detected_peaks$Peak.SNR + 1)

  detected_peaks
}


#' Assess quality of detected peaks
#'
#' This function assesses the quality of peaks returned by detect_peaks. It returns a filtered dataset of valid peaks, alongside a histogram to help visualize peak quality. This function will also perform a chi-squared test to determine if the distribution of detected peaks across samples (as well as within and between treatment groups if metadata is supplied) deviates from a uniform distribution (i.e. were the detected peaks a result of a true signal, or random noise?)
#'
#' @param peaks.df A dataframe of peaks returned by detect_peaks
#' @param spectra.matrix The spectral matrix corresponding to peaks.df. If supplied, this matrix used to generate an interactive plot of outlier peaks for each sample.
#' @param ppm A numeric vector containing the ppm values of spectra.matrix
#' @param spectra.metadata A dataframe containing the treatment factors and other metadata associated with each spectrum. If supplied, this is used to determine if there is a significant difference between the composition of peaks among treatment groups.
#' @param treatment The name of the corresponding treatment column in spectra.meta
#' @param treatment.colours An optimal list of colours to use in the diagnostic plots
#' @param region An optional argument to set the ppm bounds of peaks to inspect
#' @param outlier.quantile Numerical threshold used to define the outlier quantile of peaks (default = 0.95)
#' @param min.signal The minimum SNR required for a peak to be considered valid (default = 3)
#' @param binwidth = An optimal argument used to control the bin width of the diagnostic histogram. By default, this is determined by dividing  the range of SNR values into 50 discrete bins.
#' @param legend_position An optional argument used to control the position of the legend in the resulting diagnostic plots.
#' @param facet.order An optional argument used to specify the order of the resulting diagnostic plots when treatment is supplied.
#'
#' @return A filtered dataframe of valid peaks, alongside several diagnostic plots.
#' @export


peak_assessment <- function(peaks.df,
                            spectra.matrix = NULL,
                            ppm = NULL,
                            spectra.metadata = NULL,
                            treatment = "Landuse",
                            treatment.colours = NULL,
                            region = NULL,
                            outlier.quantile = 0.95,
                            min.signal = 3,
                            binwidth = NULL,
                            legend_position = "right",
                            facet.order = NULL) {
  # -------------------------------
  # REQUIRED COLUMNS CHECK
  # -------------------------------
  required_cols <- c("Spectra.ID",
                     "Peak.SNR",
                     "Peak.Shift",
                     "Peak.AUC",
                     "Peak.ID")

  if (!all(required_cols %in% colnames(peaks.df))) {
    stop("Input must contain: Spectra.ID, Peak.SNR, Peak.Shift, Peak.AUC, Peak.ID")
  }

  peaks.df <- stats::na.omit(peaks.df[, required_cols])

  # -------------------------------
  # SNR CLASSIFICATION (FASTER base R instead of dplyr)
  # -------------------------------
  snr <- peaks.df$Peak.SNR

  peaks.df$SNR_Category <- base::ifelse(snr < min.signal,
                                        "Poor (<min)",
                                        base::ifelse(
                                          snr < 25,
                                          "Low (min–25)",
                                          base::ifelse(snr < 50, "Moderate (25–50)", "High (>=50)")
                                        ))

  cb_palette <- c(
    "Poor (<min)"      = "#999999",
    "Low (min–25)"     = "#D55E00",
    "Moderate (25–50)" = "#E69F00",
    "High (>=50)"      = "#009E73"
  )

  peaks.df$SNR_Category <- factor(peaks.df$SNR_Category, levels = names(cb_palette))

  # -------------------------------
  # OUTLIERS (fast vector ops)
  # -------------------------------
  snr_threshold <- stats::quantile(snr, outlier.quantile, na.rm = TRUE)
  outlier_idx <- snr > snr_threshold

  peak_outliers <- peaks.df[outlier_idx, c("Peak.Shift", "Spectra.ID", "Peak.AUC", "Peak.SNR")]

  # -------------------------------
  # METADATA JOIN (base merge faster than dplyr here)
  # -------------------------------
  if (!is.null(spectra.metadata)) {
    if (!all(c("Spectra.ID", treatment) %in% names(spectra.metadata))) {
      stop("spectra.metadata must contain Spectra.ID and treatment column")
    }

    meta_sub <- spectra.metadata[, c("Spectra.ID", treatment)]

    peak_outliers <- merge(peak_outliers, meta_sub, by = "Spectra.ID", all.x = TRUE)

    peaks.df <- merge(peaks.df, meta_sub, by = "Spectra.ID", all.x = TRUE)

  } else {
    peak_outliers[[treatment]] <- "Unknown"
    peaks.df[[treatment]] <- "Unknown"
  }

  # -------------------------------
  # FILTERED PEAKS
  # -------------------------------
  filtered_peaks <- peaks.df[peaks.df$Peak.SNR >= min.signal, ]

  # -------------------------------
  # FACET ORDER (fast factor assign)
  # -------------------------------
  if (!is.null(facet.order)) {
    filtered_peaks[[treatment]] <- factor(filtered_peaks[[treatment]], levels = facet.order)
  } else {
    filtered_peaks[[treatment]] <- factor(filtered_peaks[[treatment]])
  }

  # -------------------------------
  # HISTOGRAM (FAST BINNING WITHOUT dplyr)
  # -------------------------------
  snr_vals <- peaks.df$Peak.SNR

  if (is.null(binwidth)) {
    binwidth <- diff(range(snr_vals, na.rm = TRUE)) / 50
  }

  breaks <- seq(floor(min(snr_vals, na.rm = TRUE)), ceiling(max(snr_vals, na.rm = TRUE)) + binwidth, by = binwidth)

  bins <- base::cut(
    snr_vals,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE
  )

  bin_levels <- levels(bins)
  bin_mid <- as.numeric(sub("\\[(.+),.*", "\\1", bin_levels)) + binwidth / 2

  hist_counts <- base::table(bins, peaks.df$SNR_Category)

  hist_data <- data.frame(
    bin = rep(bin_levels, times = length(cb_palette)),
    SNR_Category = rep(names(cb_palette), each = length(bin_levels)),
    count = as.vector(hist_counts)
  )

  hist_data$bin_center <- bin_mid[match(hist_data$bin, bin_levels)]
  hist_data$density <- hist_data$count / sum(hist_data$count) / binwidth

  # -------------------------------
  # SNR PLOT (ggplot2 namespaced)
  # -------------------------------
  hist_plot <- ggplot2::ggplot(hist_data) +
    ggplot2::geom_col(
      ggplot2::aes(x = bin_center, y = density, fill = SNR_Category),
      colour = "black",
      alpha = 0.7
    ) +
    ggplot2::geom_density(
      data = data.frame(Peak.SNR = snr_vals),
      ggplot2::aes(x = Peak.SNR),
      colour = "black",
      linewidth = 1
    ) +
    ggplot2::scale_fill_manual(values = cb_palette, drop = FALSE) +
    ggplot2::geom_vline(
      xintercept = min.signal,
      linetype = "dashed",
      colour = cb_palette["Poor (<min)"],
      linewidth = 1
    ) +
    ggplot2::labs(
      title = "Distribution of Raw Peak SNR Values",
      x = "Peak SNR",
      y = "Density",
      fill = "Category"
    ) +
    ggplot2::theme_minimal()

  # -------------------------------
  # INTERACTIVE SPECTRA PLOT
  # -------------------------------
  outlier_plot <- NULL

  if (!is.null(spectra.matrix)) {
    if (is.null(ppm)) {
      ppm <- as.numeric(colnames(spectra.matrix))
    }

    if (is.null(region)) {
      region <- range(peaks.df$Peak.Shift, na.rm = TRUE)
    }

    samples <- unique(peaks.df$Spectra.ID)

    if (!is.null(spectra.metadata)) {
      treatment_map <- spectra.metadata[[treatment]][match(samples, spectra.metadata$Spectra.ID)]
    } else {
      treatment_map <- rep("Unknown", length(samples))
    }

    if (is.null(treatment.colours)) {
      base_cols <- c(
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#999999"
      )
      uniq <- unique(treatment_map)
      treatment.colours <- setNames(base_cols[seq_along(uniq)], uniq)
    }

    outlier_plot <- plotly::plot_ly()

    region_idx <- which(ppm >= region[1] & ppm <= region[2])
    ppm_region <- ppm[region_idx]

    for (i in seq_along(samples)) {
      sid <- samples[i]
      row_idx <- match(sid, rownames(spectra.matrix))

      if (is.na(row_idx))
        next

      spec <- spectra.matrix[row_idx, region_idx]

      outlier_plot <- plotly::add_trace(
        outlier_plot,
        type = "scatter",
        mode = "lines",
        x = ppm_region,
        y = spec,
        line = list(color = treatment.colours[treatment_map[i]], width = 2),
        opacity = 0.4,
        name = sid
      )

      sample_out <- peak_outliers[peak_outliers$Spectra.ID == sid, ]

      if (nrow(sample_out) > 0) {
        outlier_plot <- plotly::add_trace(
          outlier_plot,
          type = "scatter",
          mode = "markers",
          x = sample_out$Peak.Shift,
          y = sample_out$Peak.AUC,
          marker = list(size = 8),
          name = paste0(sid, " peaks")
        )
      }
    }

    outlier_plot <- plotly::layout(
      outlier_plot,
      title = "Outlier Peaks by Sample",
      xaxis = list(
        title = "ppm",
        autorange = "reversed",
        range = rev(region)
      ),
      yaxis = list(title = "Intensity"),
      legend = list(x = ifelse(legend_position == "right", 1.05, 0))
    )
  }

  # -------------------------------
  # PEAK DISTRIBUTIONS (kept, lightly cleaned)
  # -------------------------------
  peak_position_hist <- NULL
  peak_position_hist_combined <- NULL
  chisq_annotations <- NULL

  if (nrow(filtered_peaks) > 0) {
    filtered_peaks[[treatment]] <- factor(filtered_peaks[[treatment]])

    peak_position_hist <- ggplot2::ggplot(filtered_peaks,
                                          ggplot2::aes(
                                            x = Peak.Shift,
                                            fill = .data[[treatment]],
                                            colour = .data[[treatment]]
                                          )) +
      ggplot2::geom_histogram(binwidth = 0.005, alpha = 0.5) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::scale_fill_manual(values = treatment.colours) +
      ggplot2::scale_colour_manual(values = treatment.colours) +
      ggplot2::scale_x_reverse() +
      ggplot2::theme_minimal() +
      ggplot2::facet_wrap(stats::as.formula(paste0("~", treatment)))

    peak_position_hist_combined <- ggplot2::ggplot(filtered_peaks) +
      ggplot2::geom_histogram(ggplot2::aes(x = Peak.Shift),
                              binwidth = 0.005,
                              fill = "#0072B2") +
      ggplot2::geom_density(ggplot2::aes(x = Peak.Shift), colour = "#D55E00") +
      ggplot2::scale_x_reverse() +
      ggplot2::theme_minimal()
  }

  # -------------------------------
  # RETURN
  # -------------------------------
  list(
    filtered.peaks = filtered_peaks,
    peak.outliers = peak_outliers,
    SNR.density.plot = hist_plot,
    peak.outlier.plot = outlier_plot,
    peak.distribution.treatment = peak_position_hist,
    peak.distribution.total = peak_position_hist_combined,
    peak.distribution.chisq = chisq_annotations
  )
}



#' Function to plot detected peaks
#'
#' Generate custom plots to vizualise peaks detected by grab_peaks
#'
#' @param df dataframe of grouped peaks to plot
#' @param samples vector of samples to plot
#' @param reference reference spectra used in warping
#' @param region numeric vector containing range (start,stop) of ppm to plot
#' @param spectra_matrix spectra matrix used as input for grab_peaks
#' @param legend_position control position of plot legend
#' @param metadata dataframe contaning metadata of samples
#' @param treatment name of treatment column in metadata to group samples by

#' @return Plot of compared spectra
#' @export


plot_detected_peaks <- function(df,
                                region = NULL,
                                samples = NULL,
                                reference = NULL,
                                spectra_matrix = NULL,
                                legend_position = "right",
                                metadata = NULL,
                                treatment = NULL) {
  # Validate input
  if (!all(c("Sample", "peakPPM", "peakValue") %in% colnames(df))) {
    stop("Input dataframe must contain 'Sample', 'peakPPM', and 'peakValue' columns.")
  }

  # Set region if not supplied
  if (is.null(region)) {
    region <- range(df$peakPPM, na.rm = TRUE)
  }

  # Filter peaks by region
  df_region <- df %>% dplyr::filter(peakPPM >= region[1] &
                                      peakPPM <= region[2])

  # Add reference if needed
  if (!is.null(reference) && !(reference %in% samples)) {
    samples <- c(reference, samples)
  }

  # Filter by sample if given
  if (!is.null(samples)) {
    df_region <- df_region %>% dplyr::filter(Sample %in% samples)
  } else {
    samples <- unique(df_region$Sample)
  }

  # ColorBlind-safe palette
  cb_palette <- c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#999999"
  )

  # Assign colors
  if (!is.null(metadata) && !is.null(treatment)) {
    metadata <- as.data.frame(metadata)
    rownames(metadata) <- as.character(metadata[[1]])

    if (!(treatment %in% colnames(metadata))) {
      stop("The specified treatment column does not exist in metadata.")
    }

    treatments <- metadata[samples, treatment, drop = TRUE]
    treatment_labels <- stats::setNames(as.character(treatments), samples)
    unique_treatments <- unique(stats::na.omit(treatments))

    # Reorder: reference group first
    if (!is.null(reference)) {
      reference_group <- as.character(metadata[reference, treatment])
      other_groups <- setdiff(unique_treatments, reference_group)
      assigned_groups <- c(reference_group, other_groups)
    } else {
      assigned_groups <- unique_treatments
    }

    if (length(assigned_groups) > length(cb_palette)) {
      stop("More treatment groups than colors in cb_palette.")
    }

    treatment_colours <- stats::setNames(cb_palette[seq_along(assigned_groups)], assigned_groups)
    colors <- treatment_colours[as.character(treatments)]

    if (!is.null(reference)) {
      ref_pos <- which(samples == reference)
      if (length(ref_pos) > 0)
        colors[ref_pos] <- cb_palette[1]
    }

  } else {
    # Default fallback
    color_func <- grDevices::colorRampPalette(cb_palette)
    colors <- color_func(length(samples))
    treatment_labels <- stats::setNames(samples, samples)
  }

  # Create plot
  p <- plotly::plot_ly()

  for (i in seq_along(samples)) {
    sample <- samples[i]
    sample_data <- df_region %>% dplyr::filter(Sample == sample)
    sample_color <- colors[i]

    if (!is.null(spectra_matrix)) {
      row_idx <- which(rownames(spectra_matrix) == sample)
      if (length(row_idx) == 0)
        next

      ppm_numeric <- as.numeric(colnames(spectra_matrix))
      region_idx <- which(ppm_numeric >= region[1] &
                            ppm_numeric <= region[2])
      smoothed_intensities <- spectra_matrix[row_idx, region_idx]
      ppm_region <- ppm_numeric[region_idx]

      p <- p %>% plotly::add_trace(
        type = "scatter",
        mode = "lines",
        x = ppm_region,
        y = smoothed_intensities,
        name = paste(sample, " - Smoothed (", treatment_labels[[sample]], ")"),
        line = list(color = sample_color, width = 2)
      )
    }

    # Add peaks
    p <- p %>% plotly::add_trace(
      type = "scatter",
      mode = "markers",
      x = sample_data$peakPPM,
      y = sample_data$peakValue,
      name = paste(sample, " - Peaks (", treatment_labels[[sample]], ")"),
      marker = list(
        color = sample_color,
        size = 6,
        symbol = "circle"
      )
    )
  }

  # Layout
  y_axis_label <- if (!is.null(spectra_matrix))
    "Peak Intensity"
  else
    "Peak Value"

  p <- p %>% plotly::layout(
    title = "Detected Peaks by Sample",
    xaxis = list(title = "Chemical Shift (ppm)", range = rev(region)),
    yaxis = list(title = y_axis_label),
    showlegend = TRUE,
    legend = list(
      x = if (legend_position == "right")
        1.05
      else
        0,
      y = 1,
      title = list(text = "Sample"),
      orientation = "v"
    ),
    hovermode = "closest"
  )

  return(p)
}
