## ──────────────────────────────────────────────────────────────
## INTEGRATION NOTE (added when wiring in the DP correspondence
## pipeline — see predict-known-peaks.R, geometry-vectors.R,
## dp-gibbs.R, etc.)
##
## This file documents summarise_reference_peaks() and
## map_reference_peaks() via roxygen comments and @export tags,
## but no function body for either is defined anywhere in this
## file or the package. As written, roxygen2::roxygenise() cannot
## attach this documentation to a function and these two names will
## NOT appear in NAMESPACE — running document() will most likely
## emit a warning (or silently skip these blocks) rather than
## export anything from this file.
##
## predict_known_peaks() (see R/predict-known-peaks.R) is the
## newly added, working replacement for the assignment stage these
## docs describe: it consumes a reference peak table
## (predicted.known.peaks, equivalent in spirit to this file's
## reference.peaks.meta.f) together with total.peaks (all detected
## candidates) and guarantees every reference peak is mapped to its
## best-matching candidate. It does not perform calibration drift
## correction, alignment-window adaptation, or the Plotly
## visualisation output described below — if that functionality is
## still needed, the actual implementations of
## summarise_reference_peaks()/map_reference_peaks() should be
## restored here (or supplied) and the relationship between the two
## assignment paths reconciled explicitly, rather than having two
## undocumented, partially-overlapping peak-assignment entry
## points in the same package.
##
## This file's existing content is left otherwise unchanged below.
## ──────────────────────────────────────────────────────────────

#' Summarise and Align Reference Peaks Across Spectra
#'
#' Detects reference peaks across multiple spectra, performs adaptive alignment using
#' internal calibration peaks, assigns chemical classes, and computes global peak
#' stability and overlap structure.
#'
#' @description
#' This function performs reference-guided peak detection and alignment for NMR-style
#' spectral matrices. It:
#' \itemize{
#'   \item Detects calibration peaks to estimate spectral drift
#'   \item Computes adaptive alignment windows per dataset
#'   \item Localises reference peaks per spectrum using constrained peak search
#'   \item Assigns peaks to chemical classes
#'   \item Builds global peak barycentres and variability estimates
#'   \item Computes peak overlap structure
#'   \item Returns per-spectrum interactive visualisations
#' }
#'
#' @param reference.peaks.meta.f data.frame. Reference peak table containing at least
#'   \code{Peak.ID} and \code{Peak.Shift}.
#'
#' @param spectra.matrix.raw numeric matrix. Raw spectra (samples x ppm bins).
#'
#' @param peak.class.map named character vector mapping \code{Peak.ID -> Class}.
#'
#' @param ppm.window.base numeric. Base ppm window for local peak search.
#'
#' @param scales integer vector. Wavelet scales (reserved for future use).
#'
#' @param smooth logical. Whether to apply moving average smoothing.
#'
#' @param k.window numeric. Multiplier for adaptive alignment window.
#'
#' @param prominence.threshold numeric. Minimum peak prominence for acceptance.
#'
#' @return A list containing:
#' \itemize{
#'   \item detections: data.frame of all detected peaks
#'   \item calibration: calibration peak estimates per spectrum
#'   \item adaptive_window: scalar alignment window
#'   \item sigma_alignment: alignment variability estimate
#'   \item global_model: barycentre statistics per peak
#'   \item transport_overlap: peak similarity matrix
#'   \item qc_summary: QC metrics per peak/class
#'   \item spectrum_plots: named list of Plotly spectra
#'   \item class_colours: mapping of class → colour
#' }
#'
#' @details
#' The method uses internal calibration peaks (A and D regions) to estimate
#' spectral drift and defines an adaptive matching window. Peaks are then
#' detected via local maximum search within this window and filtered using
#' prominence and shift constraints. A global Wasserstein-inspired barycentre
#' model is constructed to quantify cross-sample variability.
#'
#' @importFrom stats sd mad median filter
#' @importFrom dplyr bind_rows group_by summarise mutate left_join n
#' @importFrom plotly plot_ly add_lines add_markers layout
#' @export
#'
#' #' Map and Quantify Reference Peaks in Calibrated Spectra
#'
#' Performs calibration-aware peak localisation and quantification across raw
#' and standardised spectra using reference peak priors.
#'
#' @description
#' This function refines peak positions using calibration offsets and reference
#' peak priors, then quantifies peak intensities using spline-based interpolation
#' and local optimisation. It compares raw and aligned spectra to evaluate the
#' effect of spectral correction.
#'
#' @param calibrated.output output list from \code{summarise_reference_peaks()}.
#'
#' @param reference.peaks.meta.f data.frame of reference peak definitions.
#'
#' @param spectra.matrix.raw numeric matrix of raw spectra.
#'
#' @param spectra.matrix.stand numeric matrix of aligned/standardised spectra.
#'
#' @param peak.class.map named character vector mapping Peak.ID → Class.
#'
#' @param D.target numeric. Target calibration peak position (ppm).
#'
#' @param sigma.match numeric. Matching tolerance for peak localisation.
#'
#' @param smooth logical. Whether to smooth spectra prior to analysis.
#'
#' @param auc.window numeric. Integration window for AUC calculation.
#'
#' @return A list containing:
#' \itemize{
#'   \item predictions: data.frame of raw vs standardised peak metrics
#'   \item qc_summary: QC metrics summarising alignment performance
#'   \item spectrum_plots: list of paired Plotly visualisations
#'   \item class_colours: class → colour mapping
#' }
#'
#' @details
#' This method constructs spline interpolations of spectra and uses local
#' optimisation to locate peak maxima near expected chemical shifts. Peak
#' areas are computed using trapezoidal integration after baseline correction.
#' Calibration drift is corrected using internal reference peak D-shift offsets.
#'
#' @importFrom stats filter optimize sd mad
#' @importFrom dplyr bind_rows group_by summarise mutate left_join n
#' @importFrom plotly plot_ly add_lines add_markers subplot layout
#' @export
