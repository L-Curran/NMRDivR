############################################################
#' Detect peak centroid from local spectrum region
#' @export
detect_peak_centroid <- function(ppm, y) {

  baseline <- stats::median(y, na.rm = TRUE)
  idx <- which(y > baseline)

  if (length(idx) < 5) {
    return(list(ppm = NA_real_, intensity = NA_real_))
  }

  w <- y[idx] - baseline
  w[w < 0] <- 0

  if (sum(w) == 0) {
    return(list(ppm = NA_real_, intensity = NA_real_))
  }

  ppm_est <- sum(ppm[idx] * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
  intensity_est <- max(y[idx], na.rm = TRUE)

  list(ppm = ppm_est, intensity = intensity_est)
}


############################################################
#' Summarise reference peaks across spectra
#' @export
summarise_reference_peaks <- function(reference.peaks.meta.f,
                                      spectra.matrix.raw,
                                      peak.class.map,
                                      ppm.window.base = 0.04,
                                      smooth = TRUE,
                                      k.window = 3) {

  ppm <- as.numeric(colnames(spectra.matrix.raw))
  peak.ids <- unique(reference.peaks.meta.f$Peak.ID)

  mu0 <- stats::aggregate(
    Peak.Shift ~ Peak.ID,
    data = reference.peaks.meta.f,
    FUN = stats::median
  )
  names(mu0)[2] <- "mu0"

  detections <- vector("list", nrow(spectra.matrix.raw))

  for (i in seq_len(nrow(spectra.matrix.raw))) {

    y <- as.numeric(spectra.matrix.raw[i, ])

    if (smooth) {
      y <- stats::filter(y, rep(1 / 5, 5), sides = 2)
      y[is.na(y)] <- 0
    }

    spec.peaks <- list()

    for (j in seq_along(peak.ids)) {

      pid <- peak.ids[j]
      center <- mu0$mu0[mu0$Peak.ID == pid]

      idx <- which(ppm >= center - ppm.window.base &
                     ppm <= center + ppm.window.base)

      if (length(idx) < 5) next

      pk <- detect_peak_centroid(ppm[idx], y[idx])

      if (is.na(pk$ppm)) next

      spec.peaks[[length(spec.peaks) + 1]] <- data.frame(
        Spectra.ID = rownames(spectra.matrix.raw)[i],
        Peak.ID = pid,
        ppm = pk$ppm,
        mu0 = center
      )
    }

    detections[[i]] <- if (length(spec.peaks)) {
      do.call(rbind, spec.peaks)
    } else {
      NULL
    }
  }

  all <- do.call(rbind, detections)

  bary <- aggregate(
    ppm ~ Peak.ID,
    data = all,
    FUN = function(x) c(
      barycenter = stats::median(x),
      mad = stats::mad(x),
      sd = stats::sd(x),
      n = length(x)
    )
  )

  list(
    detections = all,
    global_model = bary
  )
}


############################################################
#' Map reference peaks to spectra
#' @export
map_reference_peaks <- function(calibrated.output,
                                reference.peaks.meta.f,
                                spectra.matrix.raw,
                                peak.class.map,
                                sigma.ppm = 0.01) {

  detections <- calibrated.output$detections

  mu0 <- stats::aggregate(
    Peak.Shift ~ Peak.ID,
    data = reference.peaks.meta.f,
    FUN = stats::median
  )
  names(mu0)[2] <- "mu0"

  predictions <- list()

  for (sid in unique(detections$Spectra.ID)) {

    tp <- detections[detections$Spectra.ID == sid, , drop = FALSE]

    spec.out <- list()

    for (pid in unique(tp$Peak.ID)) {

      obs <- tp[tp$Peak.ID == pid, , drop = FALSE]
      if (nrow(obs) == 0) next

      mu_ref <- mu0$mu0[mu0$Peak.ID == pid]
      ppm_i <- obs$ppm[1]

      spec.out[[length(spec.out) + 1]] <- data.frame(
        Spectra.ID = sid,
        Peak.ID = pid,
        ppm = ppm_i,
        error = ppm_i - mu_ref
      )
    }

    predictions[[sid]] <- if (length(spec.out)) {
      do.call(rbind, spec.out)
    } else {
      NULL
    }
  }

  list(
    predictions = do.call(rbind, predictions),
    qc_summary = NULL
  )
}


############################################################
#' Predict known peaks using probabilistic assignment
#' @export
predict_known_peaks <- function(total.peaks,
                                predicted.known.peaks,
                                sigma.ppm = 0.04,
                                sigma.auc = NULL,
                                use.log.auc = TRUE) {

  total.peaks <- as.data.frame(total.peaks)
  predicted.known.peaks <- as.data.frame(predicted.known.peaks)

  total.peaks$Peak.ID <- as.character(total.peaks$Peak.ID)
  total.peaks$Spectra.ID <- as.character(total.peaks$Spectra.ID)
  predicted.known.peaks$Peak.ID <- as.character(predicted.known.peaks$Peak.ID)
  predicted.known.peaks$Spectra.ID <- as.character(predicted.known.peaks$Spectra.ID)

  predicted.known.peaks <- predicted.known.peaks[
    !duplicated(predicted.known.peaks[c("Spectra.ID", "Peak.ID")]),
  ]

  if (is.null(sigma.auc)) {
    sigma.auc <- stats::sd(log1p(predicted.known.peaks$raw_auc), na.rm = TRUE) * 0.5
  }

  if (use.log.auc) {
    predicted.known.peaks$auc_used <- log1p(predicted.known.peaks$raw_auc)
    total.peaks$auc_used <- log1p(total.peaks$Peak.AUC)
  } else {
    predicted.known.peaks$auc_used <- predicted.known.peaks$raw_auc
    total.peaks$auc_used <- total.peaks$Peak.AUC
  }

  known.mu <- aggregate(
    cbind(raw_ppm, auc_used) ~ Peak.ID,
    data = predicted.known.peaks,
    FUN = mean
  )
  names(known.mu)[2:3] <- c("mu_ppm", "mu_auc")

  spectra.ids <- unique(predicted.known.peaks$Spectra.ID)

  mapped.known.list <- list()
  mapped.unknown.list <- list()
  diag.list <- list()

  for (s in seq_along(spectra.ids)) {

    sid <- spectra.ids[s]
    tp <- total.peaks[total.peaks$Spectra.ID == sid, , drop = FALSE]

    required <- predicted.known.peaks$Peak.ID[
      predicted.known.peaks$Spectra.ID == sid
    ]

    req.mu <- known.mu[known.mu$Peak.ID %in% required, , drop = FALSE]

    if (nrow(tp) == 0 || nrow(req.mu) == 0) next

    n <- nrow(tp)
    k <- nrow(req.mu)

    post.mat <- matrix(0, nrow = n, ncol = k + 1)

    for (i in seq_len(n)) {

      ppm_i <- tp$Peak.Shift[i]
      auc_i <- tp$auc_used[i]

      loglik <- -0.5 * (
        (ppm_i - req.mu$mu_ppm)^2 / sigma.ppm^2 +
          (auc_i - req.mu$mu_auc)^2 / sigma.auc^2
      )

      logpost <- c(loglik, log(1e-6))
      logpost <- logpost - max(logpost)

      post.mat[i, ] <- exp(logpost) / sum(exp(logpost))
    }

    best <- apply(post.mat[, seq_len(k), drop = FALSE], 2, which.max)

    mapped.known.list[[length(mapped.known.list) + 1]] <-
      data.frame(
        Spectra.ID = sid,
        Peak.ID = req.mu$Peak.ID,
        Mapped = TRUE
      )

    unmapped <- setdiff(seq_len(n), best)

    if (length(unmapped)) {
      mapped.unknown.list[[length(mapped.unknown.list) + 1]] <-
        tp[unmapped, , drop = FALSE]
    }

    diag.list[[length(diag.list) + 1]] <- data.frame(
      Spectra.ID = sid,
      n_peaks = n
    )
  }

  list(
    known_peaks = do.call(rbind, mapped.known.list),
    unknown_peaks = do.call(rbind, mapped.unknown.list),
    diagnostics = do.call(rbind, diag.list)
  )
}
