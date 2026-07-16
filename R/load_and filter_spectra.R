#' Load processed 1D NMR spectra into R
#'
#' This function loads preprocessed solution 1D NMR spectra from Bruker or JCAMP-DX formats
#' and returns a named list of `NMRData1D` objects, one per sample.
#'
#' @param spectra.folder Character. Path to the input folder containing spectra.
#' @param samples Character vector of sample names to load.
#' @param sample.index Character. Subfolder index within each sample directory (default "2").
#' @param procs.number Numeric. Processing file number for JCAMP-DX (optional).
#' @param blocks.number Numeric. Block number for JCAMP-DX (optional).
#' @param ntuples.number Numeric. Ntuple entry number for JCAMP-DX (optional).
#' @param format Character. Either "bruker" or "jcamp-dx".
#'
#' @return Named list of `NMRData1D` objects.
#' @export
load_NMRData1D <- function(spectra.folder,
                           samples,
                           sample.index = "2",
                           procs.number = NULL,
                           blocks.number = NULL,
                           ntuples.number = NULL,
                           format = c("bruker", "jcamp-dx")) {

  format <- match.arg(format)

  # IMPORTANT: pre-allocate full named list (keeps sample structure)
  spectra_list <- vector("list", length(samples))
  names(spectra_list) <- samples

  skipped <- character(0)

  for (sample in samples) {

    spectra_path <- file.path(spectra.folder, sample, sample.index)

    if (!base::dir.exists(spectra_path)) {
      base::message("Missing: ", sample, " -> ", spectra_path)
      skipped <- c(skipped, sample)
      spectra_list[[sample]] <- NULL   # IMPORTANT: keep slot
      next
    }

    base::message("Loading: ", sample)

    obj <- if (format == "bruker") {
      nmrdata_1d(path = spectra_path)
    } else {

      if (base::is.null(procs.number) ||
          base::is.null(blocks.number) ||
          base::is.null(ntuples.number)) {
        base::stop("JCAMP-DX requires procs.number, blocks.number, ntuples.number.")
      }

      nmrdata_1d(
        path = spectra_path,
        procs.number = procs.number,
        blocks.number = blocks.number,
        ntuples.number = ntuples.number
      )
    }

    # CRITICAL SAFETY CHECK:
    # ensure we only store valid list-like outputs
    if (!is.list(obj)) {
      base::stop("nmrdata_1d() did not return a list for sample: ", sample)
    }

    spectra_list[[sample]] <- obj
  }

  base::message(
    "Loaded ",
    sum(!vapply(spectra_list, is.null, logical(1))),
    " spectra."
  )

  if (length(skipped)) {
    base::warning("Missing samples: ", paste(skipped, collapse = ", "))
  }

  spectra_list
}
#' Ensure consistency between spectra matrix and metadata
#'
#' Filters both spectral matrix and metadata so that only shared sample IDs remain.
#'
#' @param spectra.matrix Numeric matrix (samples × ppm).
#' @param spectra.meta Data.frame containing metadata.
#' @param sample.names Column name in metadata containing sample IDs.
#'
#' @return List containing filtered spectra.matrix and spectra.meta.
#' @export
sample_sanity <- function(spectra.matrix,
                          spectra.meta,
                          sample.names = "Spectra.ID") {

  ids_matrix <- base::rownames(spectra.matrix)
  ids_meta <- unique(spectra.meta[[sample.names]])

  missing_meta <- base::setdiff(ids_matrix, ids_meta)
  if (length(missing_meta)) {
    base::warning(
      "Removing from matrix (not in metadata): ",
      paste(missing_meta, collapse = ", ")
    )
    spectra.matrix <- spectra.matrix[!ids_matrix %in% missing_meta, , drop = FALSE]
    ids_matrix <- base::rownames(spectra.matrix)
  }

  extra_meta <- base::setdiff(ids_meta, ids_matrix)
  if (length(extra_meta)) {
    base::warning(
      "Removing from metadata (not in matrix): ",
      paste(extra_meta, collapse = ", ")
    )
    spectra.meta <- spectra.meta[spectra.meta[[sample.names]] %in% ids_matrix, , drop = FALSE]
  }

  base::message(
    "Sanity check complete: ",
    length(ids_matrix), " samples retained"
  )

  list(
    spectra.matrix = spectra.matrix,
    spectra.meta = spectra.meta
  )
}
