# Cardiorespiratory coupling ----------------------------------------------
# Cardiac phase from R-peak marks, the cardiorespiratory synchrogram, and a
# high-level wrapper that reports the dominant heartbeat-per-breath locking.

#' Cardiac phase from R-peak marks
#'
#' Builds a continuous cardiac phase that advances linearly from 0 to
#' \eqn{2\pi} between consecutive R-peaks (wrapped per beat). Samples outside the
#' first/last beat are `NA`.
#'
#' @param rpeaks Integer sample indices of R-peaks (1-based).
#' @param n Length of the output phase vector (number of samples).
#' @return Numeric vector of wrapped cardiac phase in `[0, 2*pi)`, length `n`.
#' @examples
#' rpeaks <- seq(10, 1000, by = 25)
#' ph <- cardiacPhaseFromPeaks(rpeaks, 1000)
#' @export
cardiacPhaseFromPeaks <- function(rpeaks, n) {
  rpeaks <- sort(unique(as.integer(rpeaks)))
  phi <- rep(NA_real_, n)
  if (length(rpeaks) < 2L) return(phi)
  for (k in seq_len(length(rpeaks) - 1L)) {
    i0 <- rpeaks[k]; i1 <- rpeaks[k + 1L]
    if (i1 <= i0) next
    idx <- i0:(i1 - 1L)
    phi[idx] <- 2 * pi * (idx - i0) / (i1 - i0)
  }
  phi
}

#' Cardiorespiratory synchrogram
#'
#' The classic Schafer-Rosenblum synchrogram: the respiratory phase observed at
#' each heartbeat, wrapped into `m` respiratory cycles. Under n:m locking (n
#' heartbeats per m breaths) the points collapse onto `n` horizontal bands.
#'
#' @param resp_phase Instantaneous respiratory phase (length = n samples), e.g.
#'   from [instantaneousPhase()].
#' @param rpeaks Integer R-peak sample indices.
#' @param m Number of respiratory cycles to wrap into (default 1).
#' @return A data frame with `beat` (R-peak sample index) and `psi` (respiratory
#'   phase at that beat, in `[0, m)`), of class `"synchrogram"`.
#' @examples
#' t <- seq(0, 120, by = 0.01); fs <- 100
#' resp_phase <- instantaneousPhase(sin(2 * pi * 0.25 * t), fs)
#' rpeaks <- which(diff(sign(sin(2 * pi * 1.0 * t))) > 0)   # ~1 Hz beats
#' sg <- synchrogram(resp_phase, rpeaks, m = 1)
#' @export
synchrogram <- function(resp_phase, rpeaks, m = 1L) {
  n <- length(resp_phase)
  rpeaks <- sort(rpeaks[rpeaks >= 1L & rpeaks <= n])
  if (m == 1L) {
    phi <- (resp_phase[rpeaks] + 2 * pi) %% (2 * pi)   # [0, 2pi)
    psi <- phi / (2 * pi)                              # [0, 1)
  } else {
    uw <- .unwrap(resp_phase)                          # continuous phase
    psi <- (uw[rpeaks] %% (2 * pi * m)) / (2 * pi)      # [0, m)
  }
  structure(data.frame(beat = rpeaks, psi = psi),
            class = c("synchrogram", "data.frame"))
}

# minimal phase unwrap
#' @keywords internal
#' @noRd
.unwrap <- function(p) {
  if (length(p) < 2L) return(p)
  d <- diff(p)
  d <- d - 2 * pi * round(d / (2 * pi))
  c(p[1], p[1] + cumsum(d))
}

#' Analyse cardiorespiratory coupling
#'
#' High-level convenience wrapper: extracts respiratory phase (band-passed
#' Hilbert) and cardiac phase (from R-peaks), scans n:m phase locking, and
#' returns the dominant heartbeats-per-breath ratio, its locking index, and the
#' synchrogram.
#'
#' @param resp Respiration signal.
#' @param rpeaks Integer R-peak sample indices (same time base as `resp`).
#' @param sampling_rate Sampling rate in Hz.
#' @param resp_band Respiratory passband in Hz (default `c(0.1, 0.5)`).
#' @param max_ratio Largest heartbeats-per-breath ratio to consider (default 8).
#' @return An object of class `"cardioresp"`: a list with `ratio` (a
#'   `"heart:resp"` string), `n`, `m`, `lambda` (locking index), `sync_index`
#'   (1:1 phase-locking value), `synchrogram`, and the `table` of scanned ratios.
#' @examples
#' t <- seq(0, 120, by = 0.01); fs <- 100
#' resp <- sin(2 * pi * 0.25 * t)
#' card <- sin(2 * pi * 1.0 * t)                       # 4 beats per breath
#' rpeaks <- which(diff(sign(card)) > 0)
#' cr <- cardiorespiratoryCoupling(resp, rpeaks, fs)
#' cr$ratio
#' @export
cardiorespiratoryCoupling <- function(resp, rpeaks, sampling_rate,
                                      resp_band = c(0.1, 0.5), max_ratio = 8L) {
  n <- length(resp)
  phi_r <- instantaneousPhase(resp, sampling_rate, band = resp_band)
  phi_c <- cardiacPhaseFromPeaks(rpeaks, n)
  ok <- !is.na(phi_c)
  pc <- phi_c[ok]; pr <- phi_r[ok]

  # frequency ratio heart:resp = m:n, so scan orders and read off m:n
  scan <- nmPhaseSynchronization(pc, pr, orders_n = 1:2, orders_m = 1:max_ratio)
  sg <- synchrogram(phi_r, rpeaks, m = 1L)

  structure(
    list(ratio = sprintf("%d:%d", scan$m, scan$n),
         n = scan$n, m = scan$m, lambda = scan$lambda,
         sync_index = phaseSynchronizationIndex(pc, pr, 1L, 1L),
         synchrogram = sg, table = scan$table,
         sampling_rate = sampling_rate),
    class = "cardioresp")
}

#' @export
print.cardioresp <- function(x, ...) {
  cat("<Cardiorespiratory coupling>\n")
  cat(sprintf("  dominant locking : %s heartbeats:breaths (lambda = %.2f)\n",
              x$ratio, x$lambda))
  cat(sprintf("  1:1 sync index   : %.2f\n", x$sync_index))
  cat(sprintf("  heartbeats used  : %d\n", nrow(x$synchrogram)))
  invisible(x)
}

#' Plot a cardiorespiratory synchrogram
#'
#' @param x A `"synchrogram"` data frame or a `"cardioresp"` object.
#' @param sampling_rate Sampling rate, to label the x-axis in seconds (optional).
#' @param main Plot title.
#' @param ... Passed to [graphics::plot()].
#' @return `invisible(NULL)`.
#' @export
plotSynchrogram <- function(x, sampling_rate = NULL,
                            main = "Cardiorespiratory synchrogram", ...) {
  sg <- if (inherits(x, "cardioresp")) x$synchrogram else x
  stopifnot(inherits(sg, "synchrogram"))
  tt <- sg$beat
  xlab <- "beat (sample index)"
  if (!is.null(sampling_rate)) { tt <- tt / sampling_rate; xlab <- "time (s)" }
  graphics::plot(tt, sg$psi, pch = 16, cex = 0.4, col = "steelblue",
                 xlab = xlab, ylab = expression(psi),
                 ylim = c(0, max(1, max(sg$psi))), main = main, ...)
  invisible(NULL)
}
