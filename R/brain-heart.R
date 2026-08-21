# Brain-heart interaction --------------------------------------------------
# The heartbeat-evoked potential (HEP): the average brain (EEG) response
# time-locked to the R-peak, a core marker of cortical processing of cardiac
# signals. Plus a permutation test and a component-amplitude reader, and a
# lightweight directed-coupling summary that reuses the TDS lead/lag.

#' Heartbeat-evoked potential (HEP)
#'
#' Epochs a brain signal around each R-peak, baseline-corrects, and averages to
#' the heartbeat-evoked potential.
#'
#' @param eeg Numeric brain signal (one channel).
#' @param rpeaks Integer R-peak sample indices (same time base as `eeg`).
#' @param sampling_rate Sampling rate in Hz.
#' @param window Epoch window `c(pre, post)` in seconds relative to the R-peak
#'   (default `c(-0.1, 0.6)`).
#' @param baseline Baseline window `c(from, to)` in seconds for per-epoch
#'   correction; `NULL` to skip (default `c(-0.1, 0)`).
#' @return An object of class `"hep"`: a list with `time` (s, relative to the
#'   R-peak), `hep` (averaged waveform), `epochs` (epochs x time matrix), `n`
#'   (number of epochs), and `sampling_rate`.
#' @examples
#' set.seed(1); fs <- 250; n <- 30 * fs
#' rpeaks <- seq(fs, n - fs, by = round(0.9 * fs))
#' eeg <- rnorm(n)
#' hep <- heartbeatEvokedPotential(eeg, rpeaks, fs)
#' @export
heartbeatEvokedPotential <- function(eeg, rpeaks, sampling_rate,
                                     window = c(-0.1, 0.6),
                                     baseline = c(-0.1, 0)) {
  fs <- sampling_rate
  n <- length(eeg)
  pre <- round(window[1] * fs); post <- round(window[2] * fs)
  offs <- pre:post
  tvec <- offs / fs
  rp <- rpeaks[(rpeaks + pre) >= 1L & (rpeaks + post) <= n]
  if (length(rp) < 2L)
    stop("too few usable R-peaks for the requested window.", call. = FALSE)

  ep <- t(vapply(rp, function(r) eeg[r + offs], numeric(length(offs))))
  if (!is.null(baseline)) {
    bidx <- which(tvec >= baseline[1] & tvec <= baseline[2])
    if (length(bidx))
      ep <- ep - rowMeans(ep[, bidx, drop = FALSE])
  }
  structure(
    list(time = tvec, hep = colMeans(ep), epochs = ep, n = nrow(ep),
         sampling_rate = fs),
    class = "hep")
}

#' Amplitude of the HEP in a component window
#'
#' @param hep A `"hep"` object.
#' @param component `c(from, to)` seconds defining the component window
#'   (default `c(0.2, 0.4)`, a common HEP interval).
#' @param fun Summary over the window: `"mean"` (default) or `"peak"` (largest
#'   absolute value).
#' @return A single amplitude value.
#' @export
hepAmplitude <- function(hep, component = c(0.2, 0.4), fun = c("mean", "peak")) {
  stopifnot(inherits(hep, "hep"))
  fun <- match.arg(fun)
  idx <- which(hep$time >= component[1] & hep$time <= component[2])
  seg <- hep$hep[idx]
  if (fun == "mean") mean(seg) else seg[which.max(abs(seg))]
}

#' Permutation significance of the HEP component
#'
#' Tests whether the HEP component amplitude is larger than expected if the brain
#' signal were not time-locked to the heartbeat, by recomputing the component
#' amplitude for `n_perm` datasets with circularly shifted R-peaks.
#'
#' @param eeg,rpeaks,sampling_rate,window,baseline As in
#'   [heartbeatEvokedPotential()].
#' @param component Component window in seconds (see [hepAmplitude()]).
#' @param n_perm Number of permutations (default 200).
#' @return A list with `observed` (|component amplitude|), `null` (permutation
#'   amplitudes), and `p_value` (one-sided).
#' @export
hepSignificance <- function(eeg, rpeaks, sampling_rate,
                            window = c(-0.1, 0.6), baseline = c(-0.1, 0),
                            component = c(0.2, 0.4), n_perm = 200) {
  fs <- sampling_rate; n <- length(eeg)
  obs <- abs(hepAmplitude(
    heartbeatEvokedPotential(eeg, rpeaks, fs, window, baseline),
    component, "peak"))
  shifts <- round(stats::runif(n_perm, 0.05 * n, 0.95 * n))
  null <- vapply(shifts, function(s) {
    rp <- ((rpeaks + s - 1L) %% n) + 1L
    abs(hepAmplitude(
      heartbeatEvokedPotential(eeg, sort(rp), fs, window, baseline),
      component, "peak"))
  }, numeric(1))
  list(observed = obs, null = null,
       p_value = (1 + sum(null >= obs)) / (n_perm + 1))
}

#' @export
print.hep <- function(x, ...) {
  cat("<Heartbeat-evoked potential>\n")
  cat(sprintf("  epochs:     %d\n", x$n))
  cat(sprintf("  window:     %.0f to %.0f ms\n",
              1000 * x$time[1], 1000 * x$time[length(x$time)]))
  pk <- x$time[which.max(abs(x$hep))]
  cat(sprintf("  peak:       %.3f at %.0f ms\n", max(abs(x$hep)), 1000 * pk))
  invisible(x)
}

#' Plot a heartbeat-evoked potential
#'
#' @param x A `"hep"` object.
#' @param main Title.
#' @param ... Passed to [graphics::plot()].
#' @return `invisible(NULL)`.
#' @export
plotHEP <- function(x, main = "Heartbeat-evoked potential", ...) {
  stopifnot(inherits(x, "hep"))
  graphics::plot(x$time * 1000, x$hep, type = "l", lwd = 2, col = "firebrick",
                 xlab = "time from R-peak (ms)", ylab = "amplitude",
                 main = main, ...)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  graphics::abline(h = 0, col = "grey80")
  invisible(NULL)
}
