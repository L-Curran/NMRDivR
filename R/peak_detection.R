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
    include_nearbyPeaks = TRUE,
    raw_peakheight = FALSE,
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


peak_assessment <- function(peak_signals,
                            spectra_matrix = NULL,
                            ppm = NULL,
                            spectra_metadata = NULL,
                            treatment = "Landuse",
                            treatment_colours = NULL,
                            region = NULL,
                            outlier_quantile = 0.95,
                            min_signal = 3,
                            binwidth = NULL,
                            legend_position = "right",
                            facet_order = NULL) {
  # --- Check required columns ---
  required_cols <- c("Spectra.ID",
                     "Peak.SNR",
                     "Peak.Shift",
                     "Peak.AUC",
                     "Peak.ID")
  if (!all(required_cols %in% colnames(peak_signals))) {
    stop(
      "Input dataframe must contain columns: Spectra.ID, Peak.SNR, Peak.Shift, Peak.AUC and Peak.ID"
    )
  }

  peak_signals <- stats::na.omit(peak_signals[, required_cols])

  # --- Assign SNR categories ---
  cb_palette <- c(
    "Poor (<min)"       = "#999999",
    "Low (min–25)"      = "#D55E00",
    "Moderate (25–50)"  = "#E69F00",
    "High (>=50)"       = "#009E73"
  )

  peak_signals <- dplyr::mutate(
    peak_signals,
    SNR_Category = dplyr::case_when(
      Peak.SNR < min_signal ~ "Poor (<min)",
      Peak.SNR >= min_signal & Peak.SNR < 25 ~ "Low (min–25)",
      Peak.SNR >= 25 & Peak.SNR < 50 ~ "Moderate (25–50)",
      Peak.SNR >= 50 ~ "High (>=50)"
    ),
    SNR_Category = factor(SNR_Category, levels = names(cb_palette))
  )

  # --- Identify outliers ---
  snr_threshold <- stats::quantile(peak_signals$Peak.SNR, outlier_quantile)
  outlier_peaks <- dplyr::filter(peak_signals, Peak.SNR > snr_threshold)
  peak_outliers <- dplyr::select(outlier_peaks, Peak.Shift, Spectra.ID, Peak.AUC, Peak.SNR)

  # --- Metadata join ---
  if (!is.null(spectra_metadata)) {
    if (!all(c("Spectra.ID", treatment) %in% colnames(spectra_metadata))) {
      stop(
        paste0(
          "spectra_metadata dataframe must contain columns: Spectra.ID, ",
          treatment
        )
      )
    }

    peak_outliers <- dplyr::left_join(
      peak_outliers,
      dplyr::select(spectra_metadata, Spectra.ID, dplyr::all_of(treatment)),
      by = "Spectra.ID",
      relationship = "many-to-many"
    )

    peak_signals <- dplyr::left_join(
      peak_signals,
      dplyr::select(spectra_metadata, Spectra.ID, dplyr::all_of(treatment)),
      by = "Spectra.ID",
      relationship = "many-to-many"
    )
  } else {
    peak_outliers[[treatment]] <- "Unknown"
    peak_signals[[treatment]] <- "Unknown"
  }

  filtered_peaks <- dplyr::filter(peak_signals, Peak.SNR >= min_signal)

  # --- Control facet order ---
  if (!is.null(facet_order)) {
    filtered_peaks[[treatment]] <- factor(filtered_peaks[[treatment]], levels = facet_order)
  } else {
    filtered_peaks[[treatment]] <- factor(filtered_peaks[[treatment]])
  }

  # --- Histogram of SNR ---
  if (is.null(binwidth)) {
    binwidth <- (
      max(peak_signals$Peak.SNR, na.rm = TRUE) - min(peak_signals$Peak.SNR, na.rm = TRUE)
    ) / 50
  }

  hist_data <- dplyr::mutate(peak_signals,
                             bin = cut(
                               Peak.SNR,
                               breaks = seq(floor(min(Peak.SNR, na.rm = TRUE)), ceiling(max(Peak.SNR, na.rm = TRUE)) + binwidth, by = binwidth),
                               include.lowest = TRUE,
                               right = FALSE
                             )) %>%
    dplyr::group_by(bin, SNR_Category) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      bin_center = as.numeric(sub("\\[(.+),.*", "\\1", bin)) + binwidth / 2,
      density = count / sum(count) / binwidth
    )

  hist_plot <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = hist_data,
      ggplot2::aes(x = bin_center, y = density, fill = SNR_Category),
      color = "black",
      alpha = 0.7
    ) +
    ggplot2::geom_density(
      data = peak_signals,
      ggplot2::aes(x = Peak.SNR),
      color = "black",
      size = 1.2
    ) +
    ggplot2::scale_fill_manual(values = cb_palette, drop = FALSE) +
    ggplot2::geom_vline(
      xintercept = min_signal,
      linetype = "dashed",
      color = cb_palette["Poor (<min)"],
      size = 1
    ) +
    ggplot2::annotate(
      "text",
      x = min_signal,
      y = Inf,
      label = paste0("Min SNR = ", min_signal),
      vjust = -0.5,
      hjust = -0.1,
      color = cb_palette["Poor (<min)"]
    ) +
    ggplot2::labs(
      title = "Distribution of Raw Peak SNR Values",
      x = "Raw Peak SNR",
      y = "Density",
      fill = "SNR Category"
    ) +
    ggplot2::theme_minimal()

  # --- Interactive Plotly plot ---
  outlier_plot <- NULL
  if (!is.null(spectra_matrix)) {
    if (is.null(ppm))
      ppm <- as.numeric(colnames(spectra_matrix))
    if (is.null(region))
      region <- range(peak_signals$Peak.Shift, na.rm = TRUE)

    samples <- unique(peak_signals$Spectra.ID)
    if (!is.null(spectra_metadata)) {
      treatments <- setNames(as.character(spectra_metadata[[treatment]][match(samples, spectra_metadata$Spectra.ID)]), samples)
      if (!is.null(treatment_colours)) {
        colors <- treatment_colours[treatments]
      } else {
        unique_treatments <- unique(treatments)
        cb_plotly <- c(
          "#E69F00",
          "#56B4E9",
          "#009E73",
          "#F0E442",
          "#0072B2",
          "#D55E00",
          "#CC79A7",
          "#999999"
        )
        if (length(unique_treatments) > length(cb_plotly))
          stop("More treatment groups than colors")
        treatment_colours <- setNames(cb_plotly[seq_along(unique_treatments)], unique_treatments)
        colors <- treatment_colours[treatments]
      }
    } else {
      colors <- rep("#0072B2", length(samples))
    }

    outlier_plot <- plotly::plot_ly()
    for (i in seq_along(samples)) {
      sample <- samples[i]
      sample_color <- colors[i]
      row_idx <- which(rownames(spectra_matrix) == sample)
      if (length(row_idx) > 0) {
        region_idx <- which(ppm >= region[1] & ppm <= region[2])
        sample_spectrum <- spectra_matrix[row_idx, region_idx]
        ppm_region <- ppm[region_idx]

        outlier_plot <- outlier_plot %>% plotly::add_trace(
          type = 'scatter',
          mode = 'lines',
          x = ppm_region,
          y = sample_spectrum,
          line = list(color = sample_color, width = 2),
          opacity = 0.4,
          name = paste(sample, "(spectrum)")
        )
      }

      sample_outliers <- dplyr::filter(peak_outliers, Spectra.ID == sample)
      if (nrow(sample_outliers) > 0) {
        outlier_plot <- outlier_plot %>% plotly::add_trace(
          type = 'scatter',
          mode = 'markers',
          x = sample_outliers$Peak.Shift,
          y = sample_outliers$Peak.AUC,
          marker = list(
            color = sample_color,
            size = 8,
            symbol = "circle"
          ),
          name = paste(sample, "(outliers)")
        )
      }
    }

    outlier_plot <- outlier_plot %>% plotly::layout(
      title = "Outlier Peaks by Sample and Treatment",
      xaxis = list(
        title = "ppm",
        autorange = "reversed",
        range = rev(region)
      ),
      yaxis = list(title = "Peak Intensity"),
      legend = list(
        x = ifelse(legend_position == "right", 1.05, 0),
        y = 1,
        title = list(text = "<b>Sample</b>")
      ),
      hovermode = "closest"
    )
  }

  # --- Chi-squared Peak Distribution ---
  peak_position_hist <- NULL
  peak_position_hist_combined <- NULL
  chisq_annotations <- NULL
  if (nrow(filtered_peaks) > 0) {
    treatments <- levels(filtered_peaks[[treatment]])
    if (is.null(treatment_colours)) {
      cb_plotly <- c(
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#999999"
      )
      treatment_colours <- setNames(cb_plotly[seq_along(treatments)], treatments)
    }

    # Chi-squared test per treatment
    chisq_annotations <- data.frame(
      Treatment = treatments,
      label = sapply(treatments, function(trt) {
        trt_peaks <- as.numeric(filtered_peaks$Peak.Shift[filtered_peaks[[treatment]] == trt])
        if (length(trt_peaks) > 0) {
          counts <- table(cut(trt_peaks, breaks = 50))
          if (sum(counts) > 0) {
            x <- stats::chisq.test(counts)
            paste0(
              "Chi² = ",
              round(x$statistic, 2),
              ", p = ",
              format.pval(
                x$p.value,
                digits = 3,
                eps = 1e-10
              )
            )
          } else
            "No data"
        } else
          "No data"
      }),
      stringsAsFactors = FALSE
    )

    annotation_df <- dplyr::distinct(filtered_peaks, !!rlang::sym(treatment)) %>%
      dplyr::left_join(chisq_annotations, by = setNames("Treatment", treatment))

    peak_position_hist <- ggplot2::ggplot(
      filtered_peaks,
      ggplot2::aes(
        x = Peak.Shift,
        fill = !!rlang::sym(treatment),
        color = !!rlang::sym(treatment)
      )
    ) +
      ggplot2::geom_histogram(binwidth = 0.005,
                              alpha = 0.5,
                              position = "identity") +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::scale_fill_manual(values = treatment_colours) +
      ggplot2::scale_color_manual(values = treatment_colours) +
      ggplot2::scale_x_reverse() +
      ggplot2::labs(title = "Peak Distribution by Treatment", x = "ppm", y = "Count/Density") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = legend_position) +
      ggplot2::facet_wrap(as.formula(paste0("~", treatment)), scales = "free_y") +
      ggplot2::geom_text(
        data = annotation_df,
        ggplot2::aes(
          x = min(filtered_peaks$Peak.Shift),
          y = Inf,
          label = label
        ),
        inherit.aes = FALSE,
        vjust = 2,
        hjust = 1.5,
        size = 3
      )

    # Combined histogram
    combined_peaks <- as.numeric(filtered_peaks$Peak.Shift)
    if (length(combined_peaks) > 0) {
      combined_counts <- table(cut(combined_peaks, breaks = 50))
      combined_chisq <- if (sum(combined_counts) > 0)
        stats::chisq.test(combined_counts)
      else
        NULL
      combined_chisq_annotation <- data.frame(
        label = if (!is.null(combined_chisq))
          paste0(
            "Chi² = ",
            round(combined_chisq$statistic, 2),
            ", p = ",
            format.pval(
              combined_chisq$p.value,
              digits = 3,
              eps = 1e-10
            )
          )
        else
          "No data",
        x = min(filtered_peaks$Peak.Shift, na.rm = TRUE),
        y = Inf
      )

      peak_position_hist_combined <- ggplot2::ggplot(filtered_peaks, ggplot2::aes(x = Peak.Shift)) +
        ggplot2::geom_histogram(binwidth = 0.005,
                                fill = "#0072B2",
                                alpha = 0.7) +
        ggplot2::geom_density(color = "#D55E00", size = 1.2) +
        ggplot2::scale_x_reverse() +
        ggplot2::labs(title = "Combined Peak Distribution", x = "ppm", y = "Density") +
        ggplot2::theme_minimal() +
        ggplot2::geom_text(
          data = combined_chisq_annotation,
          ggplot2::aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          hjust = 1.5,
          vjust = 2,
          size = 3
        )
    }
  }

  return(
    list(
      filtered.peaks = filtered_peaks,
      peak.outliers = peak_outliers,
      SNR.density.plot = hist_plot,
      peak.outlier.plot = outlier_plot,
      peak.distribution.treatment = peak_position_hist,
      peak.distribution.total = peak_position_hist_combined,
      peak.distribution.chisq = chisq_annotations
    )
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
  if (!all(c("Spectra.ID", "Peak.PPM", "Peak.AUC") %in% colnames(df))) {
    stop("Input dataframe must contain 'Spectra.ID', 'Peak.PPM', and 'Peak.AUC' columns.")
  }

  # Set region if not supplied
  if (is.null(region)) {
    region <- range(df$Peak.PPM, na.rm = TRUE)
  }

  # Filter peaks by region
  df_region <- df %>% dplyr::filter(Peak.PPM >= region[1] &
                                      Peak.PPM <= region[2])

  # Add reference if needed
  if (!is.null(reference) && !(reference %in% samples)) {
    samples <- c(reference, samples)
  }

  # Filter by sample if given
  if (!is.null(samples)) {
    df_region <- df_region %>% dplyr::filter(Spectra.ID %in% samples)
  } else {
    samples <- unique(df_region$Spectra.ID)
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
      x = sample_data$Peak.PPM,
      y = sample_data$Peak.AUC,
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
