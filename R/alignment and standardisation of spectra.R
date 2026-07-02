#' Perform spectrum alignment using PepsNMR warping
#'
#' This function wraps `PepsNMR::Warping` with default parameters tuned for
#' complex environmental spectra (e.g., soils), which typically exhibit large
#' chemical shift variability.
#'
#' It performs spectral alignment, evaluates structural distortion using
#' Bray-Curtis distances and Procrustes analysis, and returns diagnostic plots.
#'
#' @param spectra_matrix_raw A numeric matrix of unaligned spectra.
#' Rows correspond to samples and columns correspond to ppm values.
#'
#' @param normalization.type Normalisation method passed to `PepsNMR::Warping`.
#'
#' @param reference.choice Reference selection strategy for warping.
#'
#' @param nperm Number of permutations for Procrustes significance testing.
#'
#' @param returnWarpFunc Logical; whether to return the warping function.
#'
#' @param returnReference Logical; whether to return the reference spectrum.
#'
#' @param ... Additional arguments passed directly to `PepsNMR::Warping`.
#'
#' @return A list containing:
#' \item{SpectraRaw}{Original spectral matrix.}
#' \item{SpectraWarped}{Warped spectral matrix.}
#' \item{ppm_grid}{Original ppm axis.}
#' \item{Warping}{Raw output from `PepsNMR::Warping`.}
#' \item{Procrustes}{Procrustes analysis object.}
#' \item{ProcrustesTest}{Permutation test of Procrustes fit.}
#' \item{ProcrustesPlot}{Diagnostic ggplot object.}
#'
#' @export
warp_me <- function(spectra_matrix_raw,
                    normalization.type = "mean",
                    reference.choice = "after",
                    nperm = 999,
                    returnWarpFunc = TRUE,
                    returnReference = TRUE,
                    ...) {
  args <- list(...)

  # -----------------------------
  # Palette helper (colourblind safe)
  # -----------------------------
  .get_palette <- function(n) {
    base_cb <- c("#0072B2",
                 "#D55E00",
                 "#009E73",
                 "#F0E442",
                 "#56B4E9",
                 "#E69F00",
                 "#CC79A7")

    if (n <= length(base_cb)) {
      return(base_cb[seq_len(n)])
    } else {
      warning("Requested colours exceed colourblind-safe palette. Using viridis.")
      return(viridis::viridis(n))
    }
  }

  # -----------------------------
  # PPM grid
  # -----------------------------
  ppm_grid <- colnames(spectra_matrix_raw)

  # -----------------------------
  # Warp spectra
  # -----------------------------
  warp_output <- do.call(PepsNMR::Warping, c(
    list(
      Spectrum_data = spectra_matrix_raw,
      normalization.type = normalization.type,
      reference.choice = reference.choice,
      verbose = TRUE
    ),
    args
  ))

  spectra_warped <- if (is.list(warp_output) &&
                        !is.null(warp_output$Spectrum_data)) {
    warp_output$Spectrum_data
  } else {
    warp_output
  }

  # -----------------------------
  # Distances (vectorised objects)
  # -----------------------------
  bc_raw <- vegan::vegdist(spectra_matrix_raw, method = "bray")
  bc_warped <- vegan::vegdist(spectra_warped, method = "bray")

  # -----------------------------
  # PCoA
  # -----------------------------
  pcoa_raw <- stats::cmdscale(bc_raw, k = 2)
  pcoa_warped <- stats::cmdscale(bc_warped, k = 2)

  # -----------------------------
  # Procrustes
  # -----------------------------
  proc <- vegan::procrustes(pcoa_raw, pcoa_warped, symmetric = TRUE)
  proc_test <- vegan::protest(pcoa_raw, pcoa_warped, permutations = nperm)

  # -----------------------------
  # Plot data
  # -----------------------------
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

  pal <- .get_palette(length(unique(df_long$Dataset)))

  p <- ggplot2::ggplot(df_long, ggplot2::aes(X, Y, color = Dataset)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_segment(
      data = df_proc,
      ggplot2::aes(
        x = Raw_X,
        y = Raw_Y,
        xend = Warped_X,
        yend = Warped_Y
      ),
      inherit.aes = FALSE,
      colour = "grey60",
      linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(values = setNames(pal, unique(df_long$Dataset))) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(title = "Procrustes Analysis: Raw vs Warped Spectra",
                  subtitle = paste0(
                    "R² = ",
                    round(proc_test$t0, 3),
                    ", p = ",
                    signif(proc_test$signif, 3)
                  ))

  # -----------------------------
  # Output
  # -----------------------------
  out <- list(
    SpectraRaw = spectra_matrix_raw,
    SpectraWarped = spectra_warped,
    ppm_grid = ppm_grid,
    Warping = warp_output,
    Procrustes = proc,
    ProcrustesTest = proc_test,
    ProcrustesPlot = p
  )

  class(out) <- "spectra_warped"
  out
}


#' Standardise warped spectra using internal chemical shift calibration
#'
#' This function corrects residual chemical shift variation in warped spectra
#' by aligning a known reference peak (orthophosphate) to a target ppm value.
#'
#' It also evaluates structural differences between warped and standardised
#' datasets using Bray-Curtis distances, db-RDA (capscale), and PERMANOVA.
#'
#' @param ppm Numeric vector of ppm values corresponding to columns.
#'
#' @param spectra_warped Matrix of warped spectra (samples × ppm).
#'
#' @param reference Reference spectrum name used during warping.
#'
#' @param ortho_peak_range Numeric vector (length 2) defining expected ppm range
#' for the orthophosphate peak.
#'
#' @param target_ppm Numeric value specifying target alignment position.
#'
#' @param metadata Data frame containing sample metadata.
#' Must include column `Spectra.ID`.
#'
#' @param treatment Name of grouping variable in metadata.
#'
#' @param treatment_colours Optional named vector of colours.
#'
#' @param nperm Number of permutations for PERMANOVA.
#'
#' @return A list containing:
#' \item{standardized.spectra}{Corrected spectral matrix.}
#' \item{summary.stats}{Summary statistics comparing datasets.}
#' \item{distance.matrices}{Bray-Curtis distance objects.}
#' \item{rda.models}{db-RDA models for each dataset.}
#' \item{summary.plot}{Combined diagnostic plot grid.}
#' \item{treatment.permanova}{PERMANOVA test results.}
#'
#' @export
standardise_warped_spectra <- function(ppm,
                                       spectra_warped,
                                       reference,
                                       ortho_peak_range = c(4.5, 6.5),
                                       target_ppm = 6.000,
                                       metadata,
                                       treatment,
                                       treatment_colours = NULL,
                                       nperm = 999) {
  # -----------------------------
  # Palette helper
  # -----------------------------
  .get_palette <- function(groups) {
    base_cb <- c("#0072B2",
                 "#D55E00",
                 "#009E73",
                 "#F0E442",
                 "#56B4E9",
                 "#E69F00",
                 "#CC79A7")

    n <- length(groups)

    if (n <= length(base_cb)) {
      setNames(base_cb[seq_len(n)], groups)
    } else {
      warning("Too many groups for colourblind palette. Using viridis.")
      cols <- viridis::viridis(n)
      setNames(cols, groups)
    }
  }

  # -----------------------------
  # Align metadata once
  # -----------------------------
  metadata <- metadata[!is.na(metadata[[treatment]]), , drop = FALSE]
  metadata <- metadata[metadata$Spectra.ID %in% rownames(spectra_warped), , drop = FALSE]

  common <- intersect(rownames(spectra_warped), metadata$Spectra.ID)

  spectra_warped <- spectra_warped[common, , drop = FALSE]
  metadata <- metadata[match(common, metadata$Spectra.ID), , drop = FALSE]

  # -----------------------------
  # Standardisation (vectorised apply instead of loop)
  # -----------------------------
  peak_region <- ppm >= ortho_peak_range[1] &
    ppm <= ortho_peak_range[2]

  corrected_spectra <- t(apply(spectra_warped, 1, function(spectrum) {
    peak_ppm <- ppm[peak_region][which.max(spectrum[peak_region])]
    shift <- peak_ppm - target_ppm

    stats::approx(ppm - shift, spectrum, xout = ppm, rule = 2)$y
  }))

  colnames(corrected_spectra) <- ppm
  rownames(corrected_spectra) <- rownames(spectra_warped)

  # -----------------------------
  # Distances
  # -----------------------------
  warped_dist <- vegan::vegdist(spectra_warped, method = "bray")
  standardized_dist <- vegan::vegdist(corrected_spectra, method = "bray")

  # -----------------------------
  # Summary stats
  # -----------------------------
  summary_stats <- data.frame(
    Dataset = c("Warped", "Standardized"),
    Mean = c(mean(warped_dist), mean(standardized_dist)),
    Median = c(
      stats::median(warped_dist),
      stats::median(standardized_dist)
    ),
    Variance = c(stats::var(warped_dist), stats::var(standardized_dist))
  )

  print(summary_stats)

  # -----------------------------
  # Distance plot data
  # -----------------------------
  all_distances <- data.frame(
    Distance = c(as.vector(warped_dist), as.vector(standardized_dist)),
    Dataset = factor(rep(
      c("Warped", "Standardized"), each = length(warped_dist)
    ), levels = c("Warped", "Standardized"))
  )

  # -----------------------------
  # db-RDA models
  # -----------------------------
  dist_list <- list(Warped = warped_dist, Standardized = standardized_dist)

  rda_models <- lapply(dist_list, function(dist_matrix) {
    vegan::capscale(stats::as.formula(paste("dist_matrix ~", treatment)), data = metadata)
  })

  # -----------------------------
  # RDA plots
  # -----------------------------
  groups <- unique(metadata[[treatment]])
  pal <- .get_palette(groups)

  rda_plots <- lapply(names(rda_models), function(name) {
    model <- rda_models[[name]]
    sc <- vegan::scores(model, display = "sites")
    eig <- vegan::eigenvals(model)

    eig <- if (length(eig) < 2)
      c(eig[1], 0)
    else
      eig

    rda_pct <- round(100 * eig[1:2] / sum(eig), 2)

    df <- data.frame(
      RDA1 = sc[, 1],
      RDA2 = sc[, 2],
      ID = metadata$Spectra.ID,
      Treatment = metadata[[treatment]]
    )

    ggplot2::ggplot(df, ggplot2::aes(RDA1, RDA2, color = Treatment)) +
      ggplot2::geom_point(size = 3) +
      ggrepel::geom_text_repel(ggplot2::aes(label = ID), size = 3) +
      ggplot2::stat_ellipse(
        ggplot2::aes(fill = Treatment),
        geom = "polygon",
        alpha = 0.2,
        colour = NA
      ) +
      ggplot2::scale_color_manual(values = pal) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0(name, " Spectra"),
        x = paste0("RDA1 (", rda_pct[1], "%)"),
        y = paste0("RDA2 (", rda_pct[2], "%)")
      )
  })

  # -----------------------------
  # PERMANOVA
  # -----------------------------
  treatment_permanova <- vegan::adonis2(stats::as.formula(paste("standardized_dist ~", treatment)),
                                        data = metadata,
                                        permutations = nperm)

  perm_text <- capture.output(print(treatment_permanova))

  perm_plot <- cowplot::ggdraw() +
    cowplot::draw_text(paste(perm_text, collapse = "\n"),
                       family = "mono",
                       fontface = "bold")

  # -----------------------------
  # Boxplot
  # -----------------------------
  dataset_cols <- c("Warped" = "#0072B2",
                    "Standardized" = "#009E73")

  summary_boxplot <- ggplot2::ggplot(all_distances,
                                     ggplot2::aes(Dataset, Distance, fill = Dataset)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::scale_fill_manual(values = dataset_cols) +
    ggplot2::theme_minimal()

  # -----------------------------
  # Combined plot
  # -----------------------------
  rda_combined_plot <- cowplot::plot_grid(
    plotlist = c(rda_plots, list(perm_plot, summary_boxplot)),
    ncol = 2,
    labels = c("A", "B", "C", "D")
  )

  # -----------------------------
  # Output
  # -----------------------------
  list(
    standardized.spectra = corrected_spectra,
    summary.stats = summary_stats,
    distance.matrices = list(warped = warped_dist, standardized = standardized_dist),
    rda.models = rda_models,
    summary.plot = rda_combined_plot,
    treatment.permanova = treatment_permanova
  )
}


#' Plot raw, warped, and standardised spectra
#'
#' This function generates interactive Plotly visualisations comparing
#' raw, warped, and optionally standardised spectra for each sample.
#'
#' Each output plot is sample-specific and overlays available spectral states
#' to facilitate qualitative assessment of alignment and correction.
#'
#' @param spectra.matrix.raw Unaligned spectral matrix.
#'
#' @param spectra.matrix.warped Warped spectral matrix.
#'
#' @param spectra.matrix.standardised Standardised spectral matrix.
#'
#' @param warping.reference Optional reference sample to highlight across datasets.
#'
#' @param sample.names Optional vector of sample names to subset plotting.
#'
#' @return A named list of Plotly interactive plots (one per sample).
#'
#' @export
plot_spectra <- function(spectra.matrix.raw = NULL,
                         spectra.matrix.warped = NULL,
                         spectra.matrix.standardised = NULL,
                         warping.reference = NULL,
                         sample.names = NULL) {
  # -----------------------------
  # Input check
  # -----------------------------
  if (is.null(spectra.matrix.raw) &&
      is.null(spectra.matrix.warped) &&
      is.null(spectra.matrix.standardised)) {
    stop("At least one of raw, warped, or standardised spectra matrices must be provided.")
  }

  # -----------------------------
  # Palette
  # -----------------------------
  get_palette <- function(n) {
    cb <- c("#0072B2",
            "#D55E00",
            "#009E73",
            "#F0E442",
            "#56B4E9",
            "#E69F00",
            "#CC79A7")

    if (n <= length(cb)) {
      cb[seq_len(n)]
    } else {
      warning("Palette exceeded colourblind-safe set; using viridis.")
      viridis::viridis(n)
    }
  }

  # -----------------------------
  # Convert matrix → list
  # -----------------------------
  build_list <- function(mat, name) {
    if (is.null(mat))
      return(NULL)

    lapply(seq_len(nrow(mat)), function(i) {
      list(
        name = rownames(mat)[i],
        ppm = as.numeric(colnames(mat)),
        intensity = as.numeric(mat[i, ]),
        type = name
      )
    })
  }

  raw_list    <- build_list(spectra.matrix.raw, "Raw")
  warped_list <- build_list(spectra.matrix.warped, "Warped")
  std_list    <- build_list(spectra.matrix.standardised, "standardised")

  all_spectra <- c(raw_list, warped_list, std_list)

  # -----------------------------
  # sample.names filter
  # -----------------------------
  use_multi <- !is.null(sample.names)

  if (use_multi) {
    all_spectra <- Filter(function(x)
      x$name %in% sample.names, all_spectra)
    if (length(all_spectra) == 0)
      stop("No spectra match sample.names.")
  }

  sample_names <- unique(vapply(all_spectra, `[[`, "", "name"))

  # -----------------------------
  # Colours
  # -----------------------------
  if (use_multi) {
    sample_cols <- setNames(get_palette(length(sample_names)), sample_names)
  } else {
    sample_cols <- setNames(rep("black", length(sample_names)), sample_names)
  }

  type_cols <- c(
    "Raw" = "#0072B2",
    "Warped" = "#D55E00",
    "standardised" = "#009E73"
  )

  linetypes <- c("Raw" = "dot",
                 "Warped" = "dash",
                 "standardised" = "solid")

  # -----------------------------
  # Plot list
  # -----------------------------
  plot_list <- list()

  for (s in sample_names) {
    sample_data <- Filter(function(x)
      x$name == s, all_spectra)

    p <- plotly::plot_ly()

    # -------------------------
    # sample spectra
    # -------------------------
    for (sp in sample_data) {
      p <- plotly::add_trace(
        p,
        x = sp$ppm,
        y = sp$intensity,
        type = "scatter",
        mode = "lines",
        name = sp$type,
        line = list(
          color = if (use_multi)
            sample_cols[s]
          else
            type_cols[sp$type],
          dash = linetypes[sp$type],
          width = 2
        )
      )
    }

    # -------------------------
    # REFERENCE SPECTRA (FIXED)
    # -------------------------
    if (!is.null(warping.reference)) {
      add_reference <- function(mat, label, dash_style, colour) {
        if (!is.null(mat) && warping.reference %in% rownames(mat)) {
          idx <- which(rownames(mat) == warping.reference)

          plotly::add_trace(
            p,
            x = as.numeric(colnames(mat)),
            y = as.numeric(mat[idx, ]),
            type = "scatter",
            mode = "lines",
            name = paste0("Reference (", label, ")"),
            line = list(
              color = colour,
              width = 3,
              dash = dash_style
            )
          )
        } else {
          p
        }
      }

      # raw reference
      p <- add_reference(spectra.matrix.raw, "Raw", "dot", "black")

      # warped reference
      p <- add_reference(spectra.matrix.warped, "Warped", "dash", "black")

      # standardised reference
      p <- add_reference(spectra.matrix.standardised,
                         "standardised",
                         "solid",
                         "black")
    }

    # -----------------------------
    # layout
    # -----------------------------
    p <- plotly::layout(
      p,
      title = paste0("Spectra: ", s),
      xaxis = list(title = "PPM", autorange = "reversed"),
      yaxis = list(title = "Intensity"),
      legend = list(orientation = "h")
    )

    plot_list[[s]] <- p
  }

  return(plot_list)
}
